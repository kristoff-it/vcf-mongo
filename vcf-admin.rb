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
	:force   => false,
}


OptionParser.new do |opts|
  opts.banner = "Manage VCF collectons.\n Usage: vcf-admin.rb [options] {list | rename | copy | delete | check | fix}  [collection]"
  
  opts.on("--force",  
  	"Perform destructive operations without asking for confirmation (delete, fix).") do |force|
    options[:force] = force
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


end.parse!

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

# Execute the command:
if ARGV.length == 0
	abort("Must be called with a command. Valid commands are: \n list, rename, delete, check, fix ")
end

case ARGV[0]
when 'list'
	if ARGV.length == 1 
		# List all collections:
		collections = list_collections(dbconn)
		puts "Collection list:\n"
		puts collections.join("\n")
		puts "\n"
	elsif ARGV.length == 2
		# List details about a single collection:
		begin
			metadata = collection_details(dbconn, ARGV[1])
			puts "Collection #{metadata['_id']}:"
			puts "\n"
			puts "Is consistent:    #{metadata['consistent']}"
			puts "Creation date:    #{metadata['created']}"
			puts "VCF files:"
			metadata['vcfs'].each_with_index do |file, index|
				puts "           #{index}: #{file}"
			end
			puts "Samples:"
			metadata['samples'].each do |sample|
				puts "           #{sample['name']}: in VCF \##{sample['vcfid']}"
			end
			puts "Call the `list` command with a collection name as parameter to check its details.\n"
		rescue => ex
			abort("Unable to list the collection: #{ex.message}")
		end
	else
		abort("The `list` command has 1 optional positional argument. \nUsage: list [COLLECTION] \nWhen a collection is specified the command print its details.")
	end

when 'rename'
	if ARGV.length != 3
		abort("The `rename` command requires 2 positional arguments. \nUsage: rename COLLECTION NEWNAME")
	end
	begin
		rename_collection(dbconn, ARGV[1], ARGV[2])
		puts "Collection renamed successfully."
	rescue => ex
		abort("Unable to rename the collection: #{ex.message}")
	end

# when 'copy'
# 	if ARGV.length != 3
# 		abort("The `copy` command requires 2 positional arguments. \nUsage: copy COLLECTION NEWCOLLECTION")
# 	end
# 	collection = ARGV[1]
# 	newcopy = ARGV[2]
# 	begin
# 		copy_collection(dbconn, collection, newcopy)
# 		puts "Collection copied successfully."
# 	rescue => ex
# 		abort("Unable to copy the collection: #{ex.message}")
# 	end

when 'delete'
	if ARGV.length != 2
		abort("The `delete` command requires 1 positional argument. \nUsage: delete COLLECTION [--force]")
	end
	if not options[:force]
		puts "!! WARNING: THIS OPERATION IS NOT REVERSIBLE !!"
		puts "Please type again the name of the collection you want to delete:\n"
		test = gets.chomp
		if test != ARGV[1]
			abort("\nCollection names do not match, aborting.\n\n")
		end
	end
	begin 
		delete_collection(dbconn, ARGV[1])
		puts "Collection deleted successfully."
	rescue => ex
		abort("Unable to delete the collection: #{ex.message}")
	end

when 'check'
	if ARGV.length != 1
		abort("The `check` command requires 1 positional argument. \nUsage: check COLLECTION")
	else 
		collection_list = bad_collections(dbconn)
		puts "Listing all collections in an inconsistent state:"
		collection_list.each do |collection|
			puts "Collection: #{collection[:name]} \nReason: #{collection[:reason]}"
			puts "\n"
		end
	end

when 'fix'
	if ARGV.length != 2
		abort("The `fix` command requires 1 positional argument. \nUsage: fix COLLECTION [--force]")
	else
		case check_collection(dbconn, ARGV[1])
		when :new
			abort("This collection does not exist.\n")
		when :consistent
			puts "This collection is consistent. Nothing to do.\n"
		when :spuriousM
			puts "This collection has spurious metadata remaining."
			puts "It will now be deleted.\n"
			begin
				fix_collection(dbconn, ARGV[1], :spuriousM)
			rescue => ex
				abort("Unable to fix the collection: #{ex.message}\n")
			end
		when :INIT
			puts "This collection has not completed its import operation."
			puts "It will now be deleted.\n"
			if not options[:force]
				puts "!!           WARNING: THIS OPERATION IS NOT REVERSIBLE           !!"
				puts " --> PLEASE MAKE SURE THE IMPORT OPERATION IS NOT STILL RUNNING <--"
				puts "Then type again the name of the collection you want to delete:\n"
				test = gets.chomp
				if test != ARGV[1]
					abort("\nCollection names do not match, aborting.\n\n")
				end
			end
			begin
				fix_collection(dbconn, ARGV[1], :INIT)
			rescue => ex
				abort("Unable to fix the collection: #{ex.message}")
			end
		when :APPEND
			puts "This collection has not completed its latest append import operation."
			puts "It will now be reverted by removing the partially-imported VCF files."
			puts "(if you stop the revert operation, when called again, the script will"
			puts "resume from where it left off, in other words no progress is lost)\n"

			if not options[:force]
				puts "!!           WARNING: THIS OPERATION IS NOT REVERSIBLE           !!"
				puts " --> PLEASE MAKE SURE THE IMPORT OPERATION IS NOT STILL RUNNING <--"
				puts "Then type again the name of the collection you want to revert:\n"
				test = gets.chomp
				if test != ARGV[1]
					abort("\nCollection names do not match, aborting.\n\n")
				end
			end
			begin
				fix_collection(dbconn, ARGV[1], :APPEND)
			rescue => ex
				abort("Unable to fix the collection: #{ex.message}")
			end
		end
	end
	# If we haven't aborted:
	puts "Operation completed successfully.\n\n"

else
	abort("Invalid command. Valid commands are: \n list, rename, delete, check, fix")
end




