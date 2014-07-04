require 'java'
require 'sam-1.113.jar'
require 'tribble-1.113.jar'
require 'variant-1.113.jar'
require 'mongo'
require 'pp'

include Mongo

java_import "org.broadinstitute.variant.vcf.VCFFileReader"


class VCFFileReader # obligatory monkey patching
	def walk_next
		if self.iterator.hasNext
			return self.iterator.next
		end
		return nil
	end
end


def direct_import(db, collection, vcf_filenames)
	# checks
	coll = MongoClient.new.db(db).collection(collection)
	headers, parsers, files = init_parsers(vcf_filenames)
	
	begin
		bulk = coll.initialize_ordered_bulk_op

  		count = 0
  		start = Time.now

  		walk_together(parsers).each do |row|
  			bulk.insert({'item' => row[0][1].toString})
  			count += 1
  			if count == 100
  				bulk.execute
  				count = 0
  				puts 100/(Time.now - start)
  				start = Time.now
  			end
  		end

	rescue => ex
	  puts ex
	end

	puts "collection count: #{coll.count}"
end


def init_parsers(vcf_filenames, ignore_bad_info=false)

	files = vcf_filenames.map {|f| java.io.File.new(f)}
	parsers = files.map {|f| VCFFileReader.new(f, false)}
	headers = parsers.map {|p| p.getFileHeader()}

	return headers, parsers, files #samples
end

def walk_together(parsers)
	record_buffer = parsers.map {|p| p.walk_next}
	exhausted_parsers = record_buffer.count(nil)
	Enumerator.new do |g|
		while exhausted_parsers < parsers.length do
			selected_records = []
			selected_records_ids = []
			lowest_chrom = '~~~~~~~~~~~~~~~~~~~~~~~~~~'
	 		# ~ is the highest ascii valued printable character, lexicographical max_value somehow
	 		lowest_pos = Float::INFINITY

	 		record_buffer.each_with_index do |record, index|
	 			if record == nil
	 				return
	 			end

	 			chrom = record.getChr
	 			pos = record.getStart

	 			if (chrom == lowest_chrom and pos < lowest_pos) or chrom < lowest_chrom
					selected_records = [record]
					selected_records_ids = [index]
					lowest_chrom = chrom
					lowest_pos = pos

				elsif chrom == lowest_chrom and pos == lowest_pos
					selected_records << record
					selected_records_ids << index
				end
			end


			selected_records_ids.each do |id|
				new_record = parsers[id].walk_next
				if new_record == nil
					exhausted_parsers += 1
				end
				record_buffer[id] = new_record
			end

			g.yield selected_records_ids.zip(selected_records)
		end
	end
end
