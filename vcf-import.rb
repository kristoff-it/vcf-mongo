require 'optparse'
require 'vcfdb-core.rb'

options = {
	:address         => 'localhost',
	:port            => 2213,
	:db              => 'VCF',
	:append          => false,
	:no_progress     => false,
	:chunk_size      => 1000,
	:ignore_bad_info => false
}


OptionParser.new do |opts|
  opts.banner = "Load VCF files into a collection.\n Usage: vcf-import.rb [options] collection file1.vcf file2.vcf ..."

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
  opts.on("--append",  
  	"Add the samples to a collection that might already have items inside (is slower than adding items to a new collection).") do |append|
    options[:append] = append
  end
  opts.on("--no-progress",  
  	"Disables showing of loading percentage completion, useful to remove clutter when logging stdout.") do |hide|
    options[:no_progress] = hide
  end
  opts.on("--chunk-size SIZE", Integer, 
  	"Select how many records to insert into MongoDB at a time when doing the initial import. Higher values might improve speed at the expense of memory usage. Defaults to 1000.") do |chunk|
    options[:chunk_size] = chunk
  end
  opts.on("--ignore-bad-info", 
  	"When specified, info fields that fail to respect their field definition (for example by having a string value inside an `Integer` field) are dropped with a warning. Other INFO fields from the same record are preserved if well formed.") do |ignore|
    options[:ignore_bad_info] = ignore
  end
end.parse!

# headers, parsers, files = init_parsers(ARGV)

# walk_together(parsers).each do |x|
# 	puts 'yo:', x
# end

direct_import('VCF', 'gnegne', ARGV)

puts options