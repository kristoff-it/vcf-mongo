require 'set'
require 'optparse'

# Picard parser:
require 'sam-1.113.jar'
require 'tribble-1.113.jar'
require 'variant-1.113.jar'

# Import the parser's java class:
java_import "org.broadinstitute.variant.vcf.VCFFileReader"

# MongoDB:
require 'mongo'
include Mongo

# Utils:
require 'vcfdb-core.rb'


options = {
	:address            => 'localhost',
	:port               => 27017,
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
  #debug
  opts.on("--queues",  
    "Show queues") do |queues|
    options[:queues] = queues
  end

  # MongoDB connection options
  opts.on("--address HOSTNAME",  
  	"Host address where the MongoDB instance is running. Defaults to `localhost`") do |host|
    options[:address] = host
  end
  opts.on("--port PORT", Integer,  
  	"Port number where the MongoDB instance is listening. Defaults to `?????`.") do |port|
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
  	"Number of threads that perform the merging operations. Defaults to 2.") do |merger|
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
# P1 P2 PN     # parser threads
#  |  |  | 
# b1 b2 bN     # parser buffers
#   \ | /
#  ALIGNER     # aligner thread
#     |
#     mb       # merger buffer 
#     | 
#   / | \
#  M1 M2 MN    # merger threads
#   \ | /
#     |
#    mdbb      # MongoDB buffer
#     |
#    / \
#  MDB1 MDBN   # MongoDB import threads 

# Load parsers:
files, parsers, headers, samples = load_parsers(vcf_filenames)

# Check sample names:
if samples.flatten.length != samples.flatten.to_set.length
	abort("Some sample names are colliding, aborting.")
end

# Check database status:
if options[:db] == 'VCF'
  puts "# Defaulting to `VCF` database."
end

dbconn = MongoClient.new(options[:address], options[:port]).db(options[:db])
dbstatus = check_db(dbconn)
case dbstatus
when :bad
  abort("Database `#{options[:db]}` doesn't seem to belong to this application.")
when :empty
  puts "# Empty database, performing initialization."
  init_db(dbconn)  
end

# Check collection status:
collectionstatus = check_collection(dbconn, collection)

case collectionstatus
when :new
  if options[:append]
    puts '# Collection does not exist, switching to direct import mode.'
  end
  options[:append] = false
when :consistent
  if not options[:append]
    abort("The collection already exists but you didn't specify the `--append` flag.")
  end
else
  abort('The collection is in an inconsistent state. Use vcf-admin to check (and fix) it.')
end

# Insert or update metadata inside the DB:
if options[:append]
  oldcount = update_metadata(dbconn, collection, vcf_filenames, headers, samples)
else
  oldcount = init_metadata(dbconn, collection, vcf_filenames, headers, samples)
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
		end
    # The aligner thread uses the :done element to check for end of queue:
		parser_buffers[index] << :done
	end
end

# Launch aligner thread:
aligner_thread = Thread.new do
	VCFRecordAligner.new(parser_buffers).each {|r| merger_buffer << r}	
  # Must add a :done for each merger thread:
	options[:merger_threads].times {merger_buffer << :done}
end

# Launch merger threads:
merger_threads = []
options[:merger_threads].times do |index|
	merger_threads << Thread.new do
		while (elem = merger_buffer.pop) != :done
			mongo_buffer << merge_records(elem, samples)
		end
    # Each upload thread must see a :done for each merger thread:
		options[:mongo_threads].times { mongo_buffer << :done}
	end
end

# Launch mongo threads:
mongo_threads_done = Counter.new # threadsafe counter
# (the counter exists to prevent the case where all threads die and we don't mark it as a failure)
mongo_threads = []
options[:mongo_threads].times do 
	mongo_threads << Thread.new do
    if not options[:append]
      mongo_direct_import(dbconn.collection(collection), mongo_buffer, options)
    else
      mongo_append_import(dbconn.collection(collection), mongo_buffer, options)
    end
    mongo_threads_done.increment!
	end 
end


# Wait for all threads to exit:
if options[:queues]
  puts "Displaying queues saturation status:"
  puts "(the order is: parsed-records queue, aligned-records queue, merged-records queue)"
  while mongo_threads.any? {|t| t.alive?} # TODO: check if condition is correct
	  print "\rP: #{parser_buffers[0].length} M: #{merger_buffer.length} DB: #{mongo_buffer.length} -- "
	  $stdout.flush
  end
else
  mongo_threads.each {|t| t.join}
end

# Check if all mongo threads are ok and finalize the import operation:
if mongo_threads_done.total == options[:mongo_threads]
  if options[:append]
    puts "Done importing new data, normalizing untouched records."
    update_untouched_records(dbconn, oldcounts, samples)
  end
  flag_as_consistent(dbconn, collection)
  puts "Import operations ended correctly."
else
  puts "Threads exited without aknowledging a successful completion." 
  puts "The collection has been left in an inconsistent state."
  puts "Use vcf-admin to check (and fix) it."
end

# Bye!
puts 'Done.'
