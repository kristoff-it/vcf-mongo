require 'java'
require 'sam-1.113.jar'
require 'tribble-1.113.jar'
require 'variant-1.113.jar'
require 'mongo'
require 'pp'

include Mongo

java_import "org.broadinstitute.variant.vcf.VCFFileReader"


class VCFRecordAligner < Enumerator
	def initialize (input_buffers) 
		@input_buffers = input_buffers
	end

	def compare (chrom, pos, current_chrom, current_pos)
		if (current_pos == nil) or (chrom == current_chrom and pos < current_pos) or (chrom < current_chrom)
		 	return -1
		elsif chrom == current_chrom and pos == current_pos
			return 0
		end
		return 1
	end

	def each(&block)
		record_buffer = @input_buffers.map {|b| b.pop}
		exhausted_parsers = record_buffer.count(:done)
		total_parsers = record_buffer.length
		current_chrom = nil
		current_pos = nil


		while exhausted_parsers < total_parsers do
			selected_records = []
			selected_records_ids = []
			current_pos = nil

	 		record_buffer.each_with_index do |record, index|
	 			if record == :done
	 				return
	 			end

	 			chrom = record.getChr
	 			pos = record.getStart

		 		case self.compare(chrom, pos, current_chrom, current_pos)
		 		when -1 
		 			current_chrom = chrom
		 			current_pos = pos
		 			selected_records = [record]
		 			selected_records_ids = [index]
		 		when 0
		 			selected_records << record
		 			selected_records_ids << index
				end
				#when 1 record is of a higher position, leave it for later
			end


			selected_records_ids.each do |id|
				record_buffer[id] = @input_buffers[id].pop
				if record_buffer[id] == :done
					exhausted_parsers += 1
				end
			end
			yield selected_records_ids.zip(selected_records)
		end
	end
end


def merge_records(tuples)
	merged_record = {
		'_id' => tuples[0][1].getChr + ':' + tuples[0][1].getStart.to_s,
		'CHROM' => tuples[0][1].getChr,
		'POS' => tuples[0][1].getStart,
		'IDs' => {},
		'REF' => tuples[0][1].getReference.getBaseString,
		'QUALs' => {},
		'FILTERs' => {},
		'INFOs' => {},
		'samples' => {},
		'snp' => tuples[0][1].isSNP
	}

	tuples.each do |id, record|
		merged_record['IDs'][id.to_s] = record.getID.to_s
		merged_record['QUALs'][id.to_s] = record.getPhredScaledQual.to_f
		merged_record['FILTERs'][id.to_s] = record.getFiltersMaybeNull.map {|f| f} #lel

		infos = {}
		record.getAttributes.map {|name, attribute| infos[name] = attribute}
		merged_record['INFOs'][id.to_s] = infos

		record.getGenotypes.each {|gt| merged_record['samples'][gt.getSampleName] = gt.getAlleles.map{|a| a.getDisplayString}}
	end

	# puts merged_record
	return merged_record
end

def load_parsers(vcf_filenames)
	files = vcf_filenames.map {|f| java.io.File.new(f)} # where is my tell()??? :(
	dotless_filenames = vcf_filenames.map {|f| f.gsub("\.", "<dot>")}
	parsers = files.map {|f| VCFFileReader.new(f, false)}
	headers = parsers.map {|p| p.getFileHeader.getMetaDataInInputOrder.map {|m| m.toString}}
	samples = parsers.map {|p| p.getFileHeader.getSampleNameToOffset.sort.map {|k,v| k}}
	return dotless_filenames, parsers, headers, samples
end

def check_db(db)
	collections = db.collection_names.reject{|c| c.start_with?("system.")} 
	if collections.length == 0
		return :empty
	end
	if collections.include? '__METADATA__'
		metadata = db.collection('__METADATA__').find_one('_id' => '__METADATA__', 'application' => 'VCFDBv1')
		if metadata	
			return :ok
		end
	end
	return :bad
end

def init_db(db)
	meta = {
		"_id" => '__METADATA__',
		"created" => Time.now,
		"application" => 'VCFDBv1'
	}
	db.collection('__METADATA__').insert(meta)
end

def check_collection(db, coll_name)
	table_exists = db.collection_names.include? coll_name
	meta = db.collection('__METADATA__').find_one('_id' => coll_name)

	if not table_exists and not meta
		return :new
	end

	if table_exists and not meta
		return :spuriousT
	end

	if meta and not table_exists
		return :spuriousM
	end

	if meta['consistent']
		return :consistent
	end

	return :inconsistent
end

def init_metadata(db, coll_name, files, headers, samples)
	samples_field = {}
	samples.zip(files).each {|sl, f| sl.each {|s| samples_field[s] = f}}
	headers_field = {}
	files.zip(headers).each {|f, h| headers_field[f] = h }
	meta = {
		'_id' => coll_name,
		'created' => Time.now,
		'vcfs' => headers_field,
		'samples' => samples_field,
		'consistent' => false,
		'last_inconsistency_reason' => ['INIT']
	}
	db.collection('__METADATA__').insert(meta)
end

def update_metadata(db, coll_name, files, headers, samples)
	#TODO: add race conditions checks
	dbmeta = db.collection('__METADATA__').find_one('_id' => coll_name)
	samples_field = dbmeta['samples']
	samples.zip(files).each {|sl, f| sl.each {|s| samples_field[s] = f}}
	headers_field = dbmeta['vcfs']
	files.zip(headers).each {|f, h| headers_field[f] = h }
	meta = {
		'vcfs' => headers_field,
		'samples' => samples_field,
		'consistent' => false,
		'last_inconsistency_reason' => ['APPEND', files]
	}
	db.collection('__METADATA__').update({'_id' => coll_name}, {'$set' => meta})
end

def flag_as_consistent(db, coll_name)
	db.collection('__METADATA__').update({'_id' => coll_name}, {"$set" => {"consistent" => true}})
end
