module VCFMongo
   module Import
      require 'set'
      require 'java'
      require 'htsjdk/sam-1.113.jar'
      require 'htsjdk/tribble-1.113.jar'
      require 'htsjdk/variant-1.113.jar'

      include Admin

      java_import "org.broadinstitute.variant.vcf.VCFFileReader"

      DATAMODEL_VERSION = 0.1

      # Uniform way of accessing error messages between Java exceptions and the Ruby ones.
      class StandardError
         def getMessage
            return self.message
         end
      end

      def init_db(db)
         meta = {
            "_id" => '__METADATA__',
            "created" => Time.now,
            "application" => "VCFDB",
            "version" => DATAMODEL_VERSION
         }
         db.collection('__METADATA__').insert(meta)
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
            'last_inconsistency_reason' => ['INIT'],
            'last_edit'                 => Time.now,
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
            'last_inconsistency_reason' => ['APPEND', files],
            'last_edit'                 => Time.now
         }
         db.collection('__METADATA__').update({'_id' => coll_name}, {'$set' => meta})
         return {:old_vcfs => dbmeta['vcfs'].length, :old_samples => dbmeta['samples'].length}
      end

      def flag_as_consistent(db, coll_name)
         db.collection('__METADATA__').update({'_id' => coll_name}, {"$set" => {"consistent" => true}})
      end


      class Record
         attr_reader :CHROM, :POS, :ID, :REF, :QUAL, :FILTER, :INFO, :meta, :samples
         def initialize (parser_record)
            @CHROM   = parser_record.getChr
            @POS     = parser_record.getStart
            @ID      = parser_record.getID
            @REF     = parser_record.getReference.getBaseString
            @QUAL    = parser_record.getPhredScaledQual
            @FILTER  = parser_record.getFilters.map {|f| f}
            @INFO    = Hash[parser_record.getAttributes.map {|name, attribute| [name, attribute] }]
            # @meta    = {'snp' => parser_record.isSNP}
            @samples = parser_record.getGenotypes.map do |call|
               
               # Standard fields:
               sample = {'GT' => call.getAlleleStrings.map{|al| al}}
               if call.hasAD
                  sample['AD'] = call.getAD.map{|ad| ad}
               end
               if call.hasDP
                  sample['DP'] = call.getDP
               end
               if call.hasGQ
                  sample['GQ'] = call.getGQ
               end
               if call.hasPL
                  sample['PL'] = call.getPL.map{|pl| pl}
               end


               # Extra fields:
               call.getExtendedAttributes.each do |k, v| 
                  if v.respond_to? :each
                     sample[k] = v.map {|x| x}
                  else
                     sample[k] = v
                  end
               end

               sample
            end
         end

         def to_s
            "<#{@CHROM}:#{@POS} - #{@samples.map{|s| s['GT'].toString}}>"
         end
      end



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

                  chrom = record.CHROM
                  pos = record.POS

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
            '_id'     => tuples[0][1].CHROM + ':' + tuples[0][1].POS.to_s,
            'CHROM'   => tuples[0][1].CHROM,
            'POS'     => tuples[0][1].POS,
            'IDs'     => [],
            'REF'    =>  tuples[0][1].REF,
            'QUALs'   => [],
            'FILTERs' => [],
            'INFOs'   => [],
            'samples' => [],
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
            merged_record['IDs']     << record.ID
            merged_record['QUALs']   << record.QUAL
            merged_record['FILTERs'] << record.FILTER
            merged_record['INFOs']   << record.INFO
            merged_record['samples'].concat(record.samples)
            if record.REF != merged_record['REF']
               record.samples.each {|s| s['_RR_'] = record.REF}
            end
         end
         return merged_record
      end

      def update_untouched_records(dbconn, coll_name, oldcounts, samples)
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
         dbconn.collection(coll_name).update(filtering_condition, update_operation, {:multi => true})
      end

      def load_parsers(vcf_filenames)
         files = vcf_filenames.map {|f| java.io.File.new(f)}
         parsers = files.map {|f| VCFFileReader.new(f, false)}
         samples = parsers.map {|p| p.getFileHeader.getSampleNameToOffset.sort.map {|k,v| k}}
         headers = parsers.map do |p|
            header = {
               'fileformat' => "",
               'contigs'    => [], 
               'formats'    => [],
               'filters'    => [],
               'infos'      => [],
               'others'     => []
            }

            p.getFileHeader.getMetaDataInSortedOrder.each do |m|  
               key = m.getKey
               case key
               when 'fileformat'
                  header['fileformat'] = m.getValue
               when 'FORMAT'
                  header['formats'] << {
                     'ID' => m.getID, 
                     'line' => m.toStringEncoding
                  }
               when 'INFO'
                  header['infos'] << {
                     'ID' => m.getID, 
                     'line' => m.toStringEncoding
                  }
               when 'FILTER'
                  header['filters'] << {
                     'ID' => m.getID, 
                     'line' => m.toStringEncoding
                  }
               when 'contig'
                  header['contigs'] << {
                     'ID' => m.getID, 
                     'line' => m.toStringEncoding
                  }
               else
                  header['others'] << {'key' => key, 'value' => m.getValue}
               end
            end

            header
         end

         return files, parsers, headers, samples
      end


      def mongo_direct_import(collection, queue, options, total_counter)
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
               total_counter.add(count)
                  count = 0
                  end
         end
         if count != 0
            bulk.execute
         end
      end

      def mongo_append_import(collection, queue, options, total_counter)
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
                     # 'meta'    => elem['meta']
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
                  total_counter.add(count)
                  count = 0
            end
         end
         if count != 0
            bulk.execute
            total_counter.add(count)
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

         def add (value)
            @mutex.synchronize { @total += value }
         end
      end
   end
end
