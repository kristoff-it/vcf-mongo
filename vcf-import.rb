require 'set'
require 'optparse'

require 'vcfdb-core.rb'


require 'sam-1.113.jar'
require 'tribble-1.113.jar'
require 'variant-1.113.jar'
require 'mongo'
require 'pp'

include Mongo

java_import "org.broadinstitute.variant.vcf.VCFFileReader"


options = {
	:address            => 'localhost',
	:port               => 2213,
	:db                 => 'VCF',
	:append             => false,
	:no_progress        => false,
	:mongo_chunk_size   => 500,
	:ignore_bad_info    => false,
	:merger_threads     => 3,
	:mongo_threads      => 2,
	:parser_buffer_size => 1000,
	:merger_buffer_size => 1000,
	:mongo_buffer_size  => 1000
}


OptionParser.new do |opts|
  opts.banner = "Load VCF files into a collection.\n Usage: vcf-import.rb [options] collection file1.vcf file2.vcf ..."

  # MongoDB connection options
  opts.on("--address HOSTNAME",  
  	"Host address where the MongoDB instance is running. Defaults to `localhost`") do |host|
    options[:address] = host
  end
  opts.on("--port PORT", Integer,  
  	"Port number where the RethinkDB instance is listening. Defaults to `?????`.") do |port|
    options[:port] = port
  end
  opts.on("--db DATABASE",  
  	"Database name where the VCF data is stored. Defaults to `VCF`.") do |db|
    options[:db] = db
  end

  # Import options
  opts.on("--ignore-bad-info", 
  	"When specified, info fields that fail to respect their field definition (for example by having a string value inside an `Integer` field) are dropped with a warning. Other INFO fields from the same record are preserved if well formed.") do |ignore|
    options[:ignore_bad_info] = ignore
  end
  opts.on("--append",  
  	"Add the samples to a collection that might already have items inside (is slower than adding items to a new collection).") do |append|
    options[:append] = append
  end
  opts.on("--no-progress",  
  	"Disables showing of loading percentage completion, useful to remove clutter when logging stdout.") do |hide|
    options[:no_progress] = hide
  end

  # Performance tweaks
  opts.on("--mongo-chunk-size SIZE", Integer, 
  	"Select how many records to insert into MongoDB at a time when doing the initial import. Higher values might improve speed at the expense of memory usage. Defaults to 500.") do |chunk|
    if chunk < 1
    	abort('mongo-chunk-size must be a value >= 1.')
    end
    options[:mongo_chunk_size] = chunk
  end
  opts.on("--merger-threads THREADNUM", Integer, 
  	"Number of threads that operate the merging operations. Defaults to 2.") do |merger|
  	if merger < 1
  		abort("merger-threads must be a value >= 1.")
  	end
    options[:merger_threads] = merger
  end
  opts.on("--mongo-threads THREADNUM", Integer,
  	"Number of threads that perform the import operations into MongoDB. Defaults to 2.") do |mongo|
    if mongo < 1
    	abort("mongo-threads must be a value >= 1.")
    end
    options[:mongo_threads] = mongo
  end
  opts.on("--parser-buffer-size BUFSIZE", Integer,
  	"Size of each parser's buffer (one parser is istantiated for each VCF file). Defaults to 1000 (records).") do |buffer|
    if buffer < 1
    	abort("parser-buffer-size must be a value >= 1.")
    end
    options[:parser_buffer_size] = buffer
  end
  opts.on("--merger-buffer-size BUFSIZE", Integer,
  	"Size of the mergers' buffer. Defaults to 1000.") do |merger|
  	if merger < 1
  		abort("merger-buffer-size must be a value >= 1.")
  	end
    options[:merger_buffer_size] = merger
  end
  opts.on("--mongo-buffer-size BUFSIZE", Integer,
  	"Sige of the mongo importers' buffer. Defaults to 1000.") do |mongo|
    if mongo < 1
    	abort("mongo-buffer-size must be a value >= 1.")
    end
    options[:mongo_buffer_size] = mongo
  end
end.parse!

# Validate remaining input:
if ARGV.length < 2
	abort("Must be called with a collection name and at least one vcf file as positional parameters.")
end

collection = ARGV[0]
vcf_filenames = ARGV[1..-1]
puts vcf_filenames
if vcf_filenames.length != vcf_filenames.to_set.length
	abort("You're trying to import the same VCF file twice.")
end


# Processing pipeline:
#
# p1 p2 p3     # parser threads
#  |  |  | 
# b1 b2 b3     # parser buffers
#   \ | /
#  ALIGNER     # aligner thread
#     |
#     mb       # merger buffer 
#     | 
#    /|\
#  m1 m2 m3    # merger threads
#     |
#    dbb       # MongoDB buffer
#     |
#    / \
#  mdb1 mdb2   # MongoDB import threads 

# Load parsers:
files, parsers, headers, samples = load_parsers(vcf_filenames)

# Check sample names:
if samples.flatten.length != samples.flatten.to_set.length
	abort("Some sample names are colliding, aborting.")
end

# Check database status:
if options[:db] == 'VCF'
	puts "# Defaulting to `VCF` database."
	coll = MongoClient.new.db('VCF').collection(collection)
end


# Instantiate queues:
parser_buffers = parsers.length.times.map {SizedQueue.new(options[:parser_buffer_size])}
merger_buffer = SizedQueue.new(options[:merger_buffer_size])
mongo_buffer = SizedQueue.new(options[:mongo_buffer_size])

# Launch parser threads:
parser_threads = []
parsers.length.times do |index|
	parser_threads << Thread.new do
		parsers[index].each do |record|
			parser_buffers[index] << record
			# puts parser_buffers.length
		end
		parser_buffers[index] << :done
	end
end

# Launch aligner thread:
aligner_thread = Thread.new do
	VCFRecordAligner.new(parser_buffers).each {|r| merger_buffer << r}	
	options[:merger_threads].times {merger_buffer << :done}
end

# Lanch merger threads:
merger_threads = []
options[:merger_threads].times do |index|
	merger_threads << Thread.new do
		while (elem = merger_buffer.pop) != :done
			mongo_buffer << merge_records(elem)
		end
		options[:mongo_threads].times { mongo_buffer << :done}
	end
end

# Launch mongo threads:
mongo_threads = []
options[:mongo_threads].times do 
	mongo_threads << Thread.new do
		coll = MongoClient.new.db('VCF').collection(collection)
		beginning = Time.now
		count = 0
		bulk = coll.initialize_ordered_bulk_op
		start = Time.now
		while (elem = mongo_buffer.pop) != :done
		  	bulk.insert(elem)
		  	count += 1
		  	if count == options[:mongo_chunk_size]
		  		bulk.execute
		  		count = 0
		 		print options[:mongo_chunk_size]/(Time.now - start)
		  		start = Time.now
		  	end
		end
		if count != 0
		  	bulk.execute
		end
		ending = Time.now
		puts "collection count: #{coll.count} in #{ending - beginning} seconds"
	end 
end


queues = true
while queues and mongo_threads.any? {|t| t.alive?}
	print "\rP: #{parser_buffers[0].length} M: #{merger_buffer.length} DB: #{mongo_buffer.length} -- "
	$stdout.flush
end
mongo_threads[0].join

puts 'Done.'



# puts options