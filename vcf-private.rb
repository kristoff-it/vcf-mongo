require 'optparse'

# MongoDB:
require 'mongo'
include Mongo

# Utils:
require 'vcfdb-core.rb'


options = {
	:address => 'localhost',
	:port    => 27017,
	:db      => 'VCF',
	:only_snps => false,
	:apply_filters => false,
	:min_qual => 0,

}


OptionParser.new do |opts|
  opts.banner = "Manage VCF collectons.\n Usage: vcf-admin.rb [options] {list | rename | copy | delete | check | fix}  [collection]"
  
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

  # Filtering options
  opts.on("--only-SNPs",  
  	"Get only private SNPs.") do |snps|
    options[:only_snps] = snps
  end
  opts.on("--apply-filters",  
  	"Get only records that have either 'PASS' or '.' as values.") do |filters|
    options[:apply_filters] = filters
  end
  opts.on("--min-qual MINQUAL",  
  	"Minimum value require for the QUAL field for a record to be considered.") do |minqual|
    options[:min_qual] = minqual
  end  
  opts.on("--against BADGROUP", Array,  
  	"Minimum value require for the QUAL field for a record to be considered.") do |badgroup|
    options[:bad_group] = badgroup
  end  
end.parse!

# Get collection parameter
if ARGV.length == 0
	abort("Must be called with a collection name as parameter.")
end
collection = ARGV[0]

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
  abort("This is an empty database, nothing to do.")
end

# Check collection status:
collectionstatus = check_collection(dbconn, collection)

if collectionstatus == :new
    abort("Collection does not exist.")
elsif collectionstatus != :consistent
	abort('The collection is in an inconsistent state. Use vcf-admin to check (and fix) it.')
end




# Execute the command:
if ARGV.length <= 3
	# Parse sample names
	good_names = nil
	bad_names = nil
	ARGV[1..2].each do |group|
		prefix = group[0]
		case prefix
		when '+'
			if good_names != nil
				abort("Only one parameter can repesent the `+` group.")
			end

			good_names = group[1..-1].split(',')
		when '#'
			if bad_names != nil
				abort("Only one parameter can represent the `-` group.")
			end
			bad_names = group[1..-1].split(',')
		end
	end

	if good_names == nil
		abort("You didn't specify the `+` group.")
	end

	# Get for each sample its id:
	begin
		good_sample_ids = get_sample_ids(dbconn, collection, good_names)
		bad_sample_ids = []
		if bad_names != nil
			bad_sample_ids = get_sample_ids(dbconn, collection, bad_names)
		end
	rescue => ex
		abort("Error while fetching the samples from the db: #{ex.message}.")
	end

	puts "GOOD SAMPLES:#{good_names}, ids: #{good_sample_ids}"
	puts "bad SAMPLES:#{bad_names}, ids: #{bad_sample_ids}"

	# Build the filter hash based on user preferences:
	filters = build_private_filters(options, good_sample_ids, bad_sample_ids)
	res = dbconn.collection(collection).find(filters)
	tot = 0
	res.each do |x|
		puts "#{x['POS']} - #{x['REF']}"
		x['samples'].each_with_index {|s, i| puts "#{i} => #{s['GT']}"}
		puts "\n"
		tot += 1
		if tot > 3
			exit()
		end
	end


else
	abort('TODO')
end









