require 'java'
require 'sam-1.113.jar'
require 'tribble-1.113.jar'
require 'variant-1.113.jar'
require 'mongo'
require 'pp'

include Mongo

java_import "org.broadinstitute.variant.vcf.VCFFileReader"

# The VCFRecordAligner performs (as one might imagine)
# the alingment of multiple VCF files in order to be able
# to perform merging operations concurrently.
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

	def each(&block) #TODO: remove yield logic to have less thread trashing 
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


def merge_records(tuples, samples)
	merged_record = {
		'_id'     => tuples[0][1].getChr + ':' + tuples[0][1].getStart.to_s,
		'CHROM'   => tuples[0][1].getChr,
		'POS'     => tuples[0][1].getStart,
		'IDs'     => [],
		'REF'     => tuples[0][1].getReference.getBaseString,
		'QUALs'   => [],
		'FILTERs' => [],
		'INFOs'   => [],
		'samples' => [],
		'snp'     => tuples[0][1].isSNP
	}
	array_index = 0
	tuples.each do |id, record|
		while id > array_index
			merged_record['IDs']     << nil
			merged_record['QUALs']   << nil
			merged_record['FILTERs'] << nil
			merged_record['INFOs']   << nil
			merged_record['samples'].concat([nil] * samples[array_index].length)
			array_index += 1
		end
		merged_record['IDs']     << record.getID
		merged_record['QUALs']   << record.getPhredScaledQual
		merged_record['FILTERs'] << record.getFilters.map {|f| f} #lel

		infos = {}
		record.getAttributes.each {|name, attribute| infos[name] = attribute}
		merged_record['INFOs']   << infos

		# Sample columns:
		record.getGenotypes.each do |call|
			# Standard fields:
			sample = {'GT' => call.getAlleleStrings}
			if call.hasAD
				sample['AD'] = call.getAD
			end
			if call.hasDP
				sample['DP'] = call.getDP
			end
			if call.hasGQ
				sample['GQ'] = call.getGQ
			end
			if call.hasPL
				sample['PL'] = call.getPL
			end
			# Extra fields:
			call.getExtendedAttributes.each {|k, v| sample[k] = v}


			merged_record['samples'] << sample
		end

	end
	return merged_record
end

def update_untouched_records(dbconn, oldcounts, samples)
	# Get all records that are missing the new samples by checking the sample list length:
	filtering_condition = {'IDs' => {'$size' => oldcounts[:old_vcfs]}}

	# Add nil for each empty slot:
	vcfnil = [nil] * oldcounts[:old_vcfs]
	update_operation = {'$push' => {
		'IDs'     => {'$each' => vcfnil},
		'QUALs'   => {'$each' => vcfnil},
		'FILTERs' => {'$each' => vcfnil},
		'INFOs'   => {'$each' => vcfnil},
		'samples' => {'$each' => [nil] * oldcounts[:old_samples]}
		}}
	dbconn.collection(collection).update(filtering_condition, update_operation, {:multi => true})
end

def load_parsers(vcf_filenames)
	files = vcf_filenames.map {|f| java.io.File.new(f)}
	parsers = files.map {|f| VCFFileReader.new(f, false)}
	headers = parsers.map {|p| p.getFileHeader.getMetaDataInInputOrder.map {|m| m.toString}}
	samples = parsers.map {|p| p.getFileHeader.getSampleNameToOffset.sort.map {|k,v| k}}
	return files, parsers, headers, samples
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
	samples_field = []
	samples.each_with_index {|sublist, i| sublist.each {|s| samples_field << {'name' => s, 'vcfid' => i}}}

	meta = {
		'_id'                       => coll_name,
		'created'                   => Time.now,
		'vcfs'                      => files,
		'headers'                   => headers,
		'samples'                   => samples_field,
		'consistent'                => false,
		'last_inconsistency_reason' => ['INIT']
	}
	db.collection('__METADATA__').insert(meta)
	return {:old_vcfs => 0, :old_samples => 0}
end

def update_metadata(db, coll_name, files, headers, samples)
	#TODO: add race conditions checks
	dbmeta = db.collection('__METADATA__').find_one('_id' => coll_name)
	samples_field = []
	samples.each_with_index {|sublist, i| sublist.each {|s| samples_field << {'name' => s, 'vcfid' => dbmeta['vcfs'].length + i}}}
	meta = {
		'vcfs'                      => dbmeta['vcfs'] + files,
		'headers'                   => dbmeta['headers'] + headers,
		'samples'                   => dbmeta['samples'] + samples_field,
		'consistent'                => false,
		'last_inconsistency_reason' => ['APPEND', files]
	}
	db.collection('__METADATA__').update({'_id' => coll_name}, {'$set' => meta})
	return {:old_vcfs => dbmeta['vcfs'].length, :old_samples => dbmeta['samples'].length}
end

def flag_as_consistent(db, coll_name)
	db.collection('__METADATA__').update({'_id' => coll_name}, {"$set" => {"consistent" => true}})
end


def mongo_direct_import(collection, queue, options)
	bulk = collection.initialize_ordered_bulk_op
	done_symbols_found = 0
	count = 0
	while done_symbols_found < options[:merger_threads]
        elem = queue.pop
        if elem == :done
          done_symbols_found += 1
          next
        end
		bulk.insert(elem)
		count += 1
		if count == options[:mongo_chunk_size]
		  	bulk.execute
		  	count = 0
		  	end
	end
	if count != 0
	 	bulk.execute
	end
end

def mongo_append_import(collection, queue, options)
	bulk = collection.initialize_ordered_bulk_op
	done_symbols_found = 0
	count = 0
	while done_symbols_found < options[:merger_threads]
        elem = queue.pop
        if elem == :done
          done_symbols_found += 1
          next
        end
		
        filtering_condition = {'_id' => elem['_id']}
        update_operation = {
        '$setOnInsert' => {
        	'CHROM'   => elem['CHROM'],
        	'POS'     => elem['POS'],
        	'REF'     => elem['REF'],
        	'snp'     => elem['snp']
        	},

        '$push' => {
        	'IDs'     => {'$each' => elem['IDs']},
        	'QUALs'   => {'$each' => elem['QUALs']},
        	'FILTERs' => {'$each' => elem['FILTERs']},
        	'INFOs'   => {'$each' => elem['INFOs']},
        	'samples' => {'$each' => elem['samples']}
        	}
        }
        bulk.find(filtering_condition).upsert.update(update_operation)

		count += 1
		if count == options[:mongo_chunk_size]
		  	bulk.execute
		  	count = 0
		end
	end
	if count != 0
	 	bulk.execute
	end	
end



class Counter
  attr_reader :total

  def initialize
    @total = 0
    @mutex = Mutex.new
  end

  def increment!
    @mutex.synchronize { @total += 1 }
  end
end
