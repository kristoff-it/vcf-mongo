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
   :bad_group => [],

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

  # Filtering options:
  opts.on("--against CONTROLGROUP", Array,  
   "Test for privates against a specific subset of the samples in this collection.") do |badgroup|
    options[:bad_group] = badgroup
  end  
  opts.on("--only-SNPs",  
   "Consider only private SNPs.") do |snps|
    options[:only_snps] = snps
  end
  opts.on("--apply-FILTERs",  
   "Exclude records that didn't pass all filters.") do |filters|
    options[:apply_filters] = filters
  end
  opts.on("--gt2ref-GQ MINGQ",  
   "Calls with a GQ value lower than MINGQ are considered to have the reference genotype.") do |mingq|
    options[:mingq] = mingq
  end  
  opts.on("--gt2ref-DP MINDP",  
   "Calls with a DP value lower than MINDP are considered to have the reference genotype.") do |mindp|
    options[:mindp] = mindp
  end  
  opts.on("--max-gt2ref-priv MAXPRIV",  
   "Records with a number higher than MAXPRIV of low quality calls in the private group are not considered.") do |maxpriv|
    options[:max_lowq_priv] = maxpriv
  end  
  opts.on("--max-gt2ref-contr MAXCONTR",  
   "Records with a number higher than MAXCONTR of low quality calls in the contol group are not considered.") do |maxcontr|
    options[:max_lowq_contr] = maxcontr
  end  
  opts.on("--max-gt2ref-total MAXLOWQ",  
   "Records with a number higher than MAXLOWQ of low quality calls (indipendently of what group they belong to) are not considered. ") do |minqual|
    options[:max_lowq_total] = minqual
  end  

end.parse!

# Check database status:
# if options[:db] == 'VCF'
#   puts "# Defaulting to `VCF` database."
# end

dbconn = MongoClient.new(options[:address], options[:port]).db(options[:db])
dbstatus = check_db(dbconn)
case dbstatus
when :bad
  abort("Database `#{options[:db]}` doesn't seem to belong to this application.")
when :empty
  abort("This is an empty database, nothing to do.")
end




# Execute the command:
if ARGV.length < 1
   abort("Must be called with a command. Valid commands are: \n get, createindex, deleteindex, serve")
end

case ARGV[0]
when "get"
   arguments = ARGV[1..-1]
   if arguments.length != 2
      abort("The `get` command requires 2 positional arguments. \nUsage: get COLLECTION PRIVATEGROUP [--against CONTROLGROUP] \n(eg: get myHumans NA01,NA02,NA03) \nWhen no control group is pecified, privates are tested against all other samples in the collection.")
   end
   
   # Check collection status:
   collection = arguments[0]
   collectionstatus = check_collection(dbconn, arguments[0])
   puts collectionstatus
   if collectionstatus == :new
       abort("Collection does not exist.")
   elsif collectionstatus != :consistent
      abort('The collection is in an inconsistent state. Use vcf-admin to check (and fix) it.')
   end

   # Parse the group argument:
   good_names = arguments[1].split(',')

   # Get the collection metadata:
   metadata = collection_details(dbconn, collection)

   # Get IDs for all samples:
   # (if a control group wasn't specified, we assume it's comprised of all the remaining samples)
   name2id = Hash[metadata['samples'].each_with_index.map {|sample, id| [sample['name'], id]}]
   puts "name2id", name2id
   good_samples = good_names.map do |name| 
      if name2id[name] != nil
         name2id[name]
      else
         abort("Sample `#{name}` is not present in this collection.")
      end
   end

   if options[:bad_group].length > 0 
      bad_samples = options[:bad_group].map do |name|
         if name2id[name] != nil
            name2id[name]
         else
            abort("Sample `#{name}` is not present in this collection")
         end
      end
   else
      bad_samples = name2id.keep_if{|k, v| not good_names.include? k}.map{|k, v| v}
   end

   relevant_vcf_ids = Set.new
   good_samples.each{|s| relevant_vcf_ids.add(metadata['samples'][s]['vcfid'])}
   bad_samples.each{|s| relevant_vcf_ids.add(metadata['samples'][s]['vcfid'])}

   # Build filters based on user preferences:
   filters = build_private_filters(:good_samples   => good_samples, 
                           :bad_samples    => bad_samples, 
                           :relevant_vcfs  => relevant_vcf_ids,
                           :only_snps      => options[:onyl_snps],
                           :apply_filters  => options[:apply_filters],
                           :mingq          => options[:mingq],
                           :mindp          => options[:mindp],
                           :max_lowq_priv  => options[:max_lowq_priv],
                           :max_lowq_contr => options[:max_lowq_contr],
                           :max_lowq_total => options[:max_lowq_total])
                           

   # Perform the query:
   res = dbconn.collection(arguments[0]).find(filters, {:fields => {'_id' => 0, 'meta' => 0}})


   # Headers
   puts '##fileformat=VCFv4.1'
   puts '##FORMAT=<ID="RR", Count=1, Type=String, Description="Relative REF: when merging multiple VCF files the REF field might not coincide (a variant might be a SNP in a VCF and a Complex in another, for example).">'

   # # INFOs
   # selected_infos = {}
   # relevant_vcf_ids.each do |vcf|
   #  metadata['headers'][vcf]['infos'].each {|i| selected_infos[i['ID']] = i['line']}
   # end
   # selected_infos.each_value {|line| puts line}

   # # FORMATs
   # selected_formats = {}
   # relevant_vcf_ids.each do |vcf|
   #  metadata['headers'][vcf]['formats'].each {|f| selected_formats[f['ID']] = f['line']}
   # end
   # selected_formats.each_value {|line| puts line}

   # Table header:
   puts "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + good_samples.map{|id| metadata['samples'][id]['name']}.join("\t")
   
   # Lines:
   res.each do |r|
      alt, format, samples = normalize(r['samples'])
      puts "#{r['CHROM']}\t#{r['POS']}\t.\t#{r['REFs'][0]}\t#{alt.join(',')}\t#{r['QUALs'].reduce(:+)}\t.\t#{format.join(':')}\t#{samples.join("\t")}"
   end

   # Present the results as a VCF file:   
   tot = 0
   res.each do |x|
      puts "#{x['POS']} - #{x['REFs']}"
      x['samples'].each_with_index {|s, i| puts "#{i} => #{s['GT']}"}
      puts "\n"
      tot += 1
      if tot > 3
         exit()
      end
   end
when "createindex"
when "deleteindex"
when "serve"
   require 'sinatra'
   get '/:collection/:samplelist' do
      "test"
   end
else
   abort("Valid commands are: get, createindex, deleteindex, serve.")
end


