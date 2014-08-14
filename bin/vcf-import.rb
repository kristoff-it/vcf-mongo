#!/usr/bin/env jruby

require 'set'
require 'optparse'


# MongoDB:
require 'mongo'
include Mongo

# Utils:
require 'lib/vcf-mongo'
include VCFMongo::Import


options = {
   :address            => 'localhost',
   :port               => 27017,
   :db                 => 'VCF',
   :append             => false,
   :no_progress        => false,
   :mongo_chunk_size   => 500,
   :merger_threads     => 1,
   :mongo_threads      => 2,
   :parser_buffer_size => 1000,
   :merger_buffer_size => 1000,
   :mongo_buffer_size  => 1000,
   :drop_bad_records   => false,
}
# Flag to check if a thread reported an error:
fatal_errors = false


OptionParser.new do |opts|
   opts.banner = "Load VCF files into a collection.\n Usage: vcf-import.rb [options] collection file1.vcf file2.vcf ..."

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

   opts.on("--append",  
      "Add the samples to a collection that might already have items inside (is slower than adding items to a new collection).") do |append|
      options[:append] = append
   end
   opts.on("--no-progress",  
      "Disables showing of loading percentage completion, useful to remove clutter when logging stdout.") do |hide|
      options[:no_progress] = hide
   end
   opts.on("--drop-bad-records",  
      "Drop malformed records without causing the whole import operation to fail. Info about dropped records is printed to STDERR.") do |non_fatal_records|
      options[:drop_bad_records] = non_fatal_records
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
   abort("Must be called with a collection name and at least one VCF filename as positional parameters.")
end

collection = ARGV[0]
if collection.include? '__'
   abort('Double underscores are used internally and cannot be part of a collection name.')
end

vcf_filenames = ARGV[1..-1]
if vcf_filenames.length != vcf_filenames.to_set.length
   # This is to prevent possibly confusing 'sample names collision' errors.
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
when :ok
   # noop
when :bad
   abort("Database `#{options[:db]}` doesn't seem to belong to this application.")
when :empty
   puts "# Empty database, performing initialization."
   init_db(dbconn)
else
   abort dbstatus
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
   oldcounts = update_metadata(dbconn, collection, vcf_filenames, headers, samples)
else
   oldcounts = init_metadata(dbconn, collection, vcf_filenames, headers, samples)
end

# Instantiate buffer queues (and error queues):
parser_buffers = parsers.length.times.map {SizedQueue.new(options[:parser_buffer_size])}
merger_buffer = SizedQueue.new(options[:merger_buffer_size])
mongo_buffer = SizedQueue.new(options[:mongo_buffer_size])

parser_threads_errors = Queue.new
aligner_thread_error = nil # since it's a single thread a queue isn't needed
merger_threads_errors = Queue.new
mongo_threads_errors = Queue.new

# Launch parser threads:
parser_threads = []
parsers.length.times do |index|
   parser_threads << Thread.new do
      begin
         parsers[index].each do |record|
            # Record is an object that mines the HTSJDK record object (called VariantContext)
            # and basically "unravels" the lazy operations. This is not optimal :sadface:,
            # but apparently different instances of VariantContext can't be used concurrently.
            # The HTSJDK parser must have gone full singleton. Never go full singleton.
            parser_buffers[index] << Record.new(record)
         end
      rescue => ex
         parser_threads_errors << ex
         fatal_errors = true
      end
      # The aligner thread uses the :done element to check for end of queue:
      parser_buffers[index] << :done
   end
end

# Launch aligner thread:
aligner_thread = Thread.new do
   begin
      VCFRecordAligner.new(parser_buffers).each {|r| merger_buffer << r}
   rescue => ex
      aligner_thread_error = ex
      fatal_errors = true
   end
   # Must add a :done for each merger thread:
   options[:merger_threads].times {merger_buffer << :done}
end

# Launch merger threads:
merger_threads = []
options[:merger_threads].times do |index|
   merger_threads << Thread.new do
      while (elem = merger_buffer.pop) != :done
         begin
            mongo_buffer << merge_records(elem, samples)
         rescue => ex
            merger_threads_errors << [elem, ex]
            if not options[:drop_bad_records]
               fatal_errors = true
               break
            end
         end
      end
      # Each upload thread must see a :done for each merger thread:
      options[:mongo_threads].times { mongo_buffer << :done}
   end
end

# Launch mongo threads:
mongo_threads_done = Counter.new # threadsafe counter
# (the counter exists to prevent the case where all threads die and we don't mark it as a failure)
imported_records_counter = Counter.new
mongo_threads = []
options[:mongo_threads].times do |index|
   mongo_threads << Thread.new do
      begin 
         if not options[:append]
            mongo_direct_import(dbconn.collection(collection), mongo_buffer, options, imported_records_counter)
         else
            mongo_append_import(dbconn.collection(collection), mongo_buffer, options, imported_records_counter)
         end
         mongo_threads_done.increment!
      rescue => ex
         mongo_threads_errors[index] = ex
         fatal_errors = true
      end
   end 
end


# Wait for all threads to exit and meanwhile show the progress status if enabled:
if not options[:no_progress]
   puts "Displaying total imported records, import speed and queues saturation status."
   puts "(queue order is: parser | merger | mongo )"
   puts "Sorry but completion percentage is not supported at the moment: https://github.com/samtools/htsjdk/issues/63\n"
end

start_timer = Time.now
previous_line_length = 0

while mongo_threads.any? {|t| t.alive?} 
   if options[:drop_bad_records]
      while not merger_threads_errors.empty?
         item, ex = merger_threads_errors.pop
         warn "\nDropping #{item[0][1].CHROM}:#{item[0][1].POS} => #{ex.getMessage}"
      end
   end

   if fatal_errors
      break
   end

   if not options[:no_progress]
      total = imported_records_counter.total
      speed = total/(Time.now + 0.001 - start_timer)
      line =  "\rTotal: #{total} @ #{speed.to_i} records/s (#{parser_buffers[0].length}|#{merger_buffer.length}|#{mongo_buffer.length})"
      if (num_spaces = (previous_line_length - line.length)) > 0
         line += ' ' * num_spaces
      end
      previous_line_length = line.length
      print line
      $stdout.flush
   end
   sleep 0.2
end

if not options[:no_progress]
   line = "\rImported #{imported_records_counter.total} records in #{Time.now - start_timer} secods."
   if (num_spaces = (previous_line_length - line.length)) > 0
      line += ' ' * num_spaces
   end
   puts line
end

puts "\n"

# There might be some more errors that must be notified:
if options[:drop_bad_records] and not fatal_errors
   while not merger_threads_errors.empty?
      item, ex = merger_threads_errors.pop
      warn "\nDropping #{item[0][1].CHROM}:#{item[0][1].POS} => #{ex.getMessage}"
   end
end

# Check if all mongo threads are ok and finalize the import operation:
if fatal_errors
   if not parser_threads_errors.empty?
      warn "Errors were encountered while executing the initial parsing operation:"
      warn "Please note that most likely these errors are bubbling up directly from the HTSJDK parser.\n"
      while parser_threads_errors.empty?
         warn parser_threads_errors.pop.getMessage
         warn "\n"
      end
   end

   if aligner_thread_error
      warn "Errors were encountered while executing the alignment phase:"
      warn aligner_thread_error.getMessage
      warn "\n"
   end

   if not merger_threads_errors.empty?
      warn "Errors were encountered while executing the merge phase:"
      while not merger_threads_errors.empty?
         item, ex = merger_threads_errors.pop
         warn "\n #{item[0][1].CHROM}:#{item[0][1].POS} => #{ex.getMessage}"
         warn "\n"
      end
   end

   if not mongo_threads_errors.empty?
      warn "Errors were encountered while executing mongo bulk insert/update operations:"
      while not mongo_threads_errors.empty?
         err = mongo_threads_errors.pop
         warn err.getMessage
         warn "\n"
      end
   end

end

if mongo_threads_done.total == options[:mongo_threads] and not fatal_errors
   if options[:append]
      puts "Done importing new data, normalizing untouched records."
      update_untouched_records(dbconn, collection, oldcounts, samples)
   end
   flag_as_consistent(dbconn, collection)
   puts "Import operations ended correctly."
   puts "Done.\n"
else
   warn "Threads exited without aknowledging a successful completion." 
   warn "The collection has been left in an inconsistent state."
   abort "Use vcf-admin to check (and fix) it.\n"
end
