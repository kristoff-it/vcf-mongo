module VCFMongo

   # The Admin module contains all administrative operations, 
   # such as listing all present collections and checking for various statuses.
   module Admin
      
      # Takes in input a Mongo::DB instance and checks if the database it references
      # belongs to this application or not. 
      # Return values are either the symbols :empty, :ok, :bad or a string that says
      # that there is a version mismatch between the script and the datamodel version
      # upon which the database was created.
      def check_db(db)
         collections = db.collection_names.reject{|c| c.start_with?("system.")} 
         if collections.length == 0
            return :empty
         end
         if collections.include? '__METADATA__'
            metadata = db.collection('__METADATA__').find_one('_id' => '__METADATA__', 'application' => 'VCFDB')
            if metadata 
               if metadata['version'] == DATAMODEL_VERSION
                  return :ok
               else
                  return "Version mismatch: script is version #{DATAMODEL_VERSION} while DB is version #{metadata['version']}."
               end
            end
         end
         return :bad
      end

      # Takes in imput a Mongo::DB instance and a String representing a collection name.
      # Checks the consistency status of the collection. Return values are the symbols
      # :new         (collection does not exist),
      # :consistent  (ok)
      # :spuriousM   (there is a metadata document but the collection does not exist),
      # :INIT        (there is an initial import operation pending)
      # :APPEND      (there is an append import operation pending)
      def check_collection(db, coll_name)
         table_exists = db.collection_names.include? coll_name
         meta = db.collection('__METADATA__').find_one('_id' => coll_name)

         if not table_exists and not meta
            return :new
         end

         if meta and not table_exists
            return :spuriousM
         end

         if meta['consistent']
            return :consistent
         end
         if meta['last_inconsistency_reason'][0] == 'INIT'
            return :INIT
         end

         return :APPEND
      end

      # Takes a Mongo::DB instance in input and returns a list of all present 
      # collections (as present in the '__METADATA__' meta-collection)
      def list_collections(db)
         return db.collection('__METADATA__').find({'_id' => {'$ne'=> '__METADATA__'}}, {:fields => '_id'}).map {|c| c['_id']}
      end

      # Takes a Mongo::DB instance and a String representing a collection name in input.
      # Returns the related metadata document containing informations about creation date,
      # imported VCF files, samples, and cosistency status.
      def collection_details(db, coll_name)
         if (result = db.collection('__METADATA__').find_one('_id' => coll_name)) == nil
            raise "collection does not exist"
         end
         return result
      end

      # Takes a Mongo::DB instance and two Strings repesenting collection names.
      # Renames an existing collection.
      def rename_collection(db, coll_name, new_name)
         collections = db.collection_names

         if new_name.include?('__')
            raise 'double underscores are resereved for internal usage'
         end

         if not collections.include? coll_name
            raise "collection does not exist"
         end
         if collections.include? new_name
            raise "target collection already exists"
         end

         db.collection(coll_name).rename(new_name)
         db.collection_names.keep_if{|c| c.start_with? coll_name + '__'}.each do |related_collection|
            db.collection(related_collection).rename(new_name + '__' + related_collection.split('__', 2)[1])
         end
         metadata = db.collection('__METADATA__').find_one('_id' => coll_name)
         metadata['_id'] = new_name
         db.collection('__METADATA__').insert(metadata)
         db.collection('__METADATA__').remove('_id' => coll_name)
      end


      # Takes a Mong::DB instance and a String representing a collection name.
      # Deletes an existing collection.
      def delete_collection(db, coll_name)
         if db.collection('__METADATA__').find_one('_id' => coll_name) == nil
            raise "collection does not exist"
         end

         db.collection(coll_name).drop
         db.collection_names.keep_if {|c| c.start_with? coll_name+'__'}.each {|related_collection| db.collection(related_collection).drop}

         db.collection('__METADATA__').remove('_id' => coll_name)
      end



      # Takes a Mongo::DB instance and returns a {:name, :reson} object representing
      # the inconsitency reason of every inconsistent collection.
      def bad_collections(db)
         return db.collection('__METADATA__').find('_id' => {'$ne' => '__METADATA__'}, 'consistent' => false).map do |c|
            if c['last_inconsistency_reason'][0] == 'INIT'
               {:name => c['_id'], :reason => 'Initial import operation non complete.'}
            else
               {:name => c['_id'], :reason => "Append import non complete (#{c['last_inconsistency_reason'][1].join(', ')})"}
            end
         end
      end

      # Takes a Mongo::DB instance, a String representing a collection name 
      # and a symbol between [:spuriousM, :INIT, :APPEND] in input.
      # Performs the appropriate fixing operation.
      def fix_collection(db, coll_name, fixtype)
         case fixtype
         when :spuriousM
            db.collection('__METADATA__').remove('_id' => coll_name)
         when :INIT
            db.collection(coll_name).drop
            db.collection('__METADATA__').remove('_id' => coll_name)
         when :APPEND
            metadata = db.collection('__METADATA__').find_one('_id' => coll_name)
            bad_vcfs = metadata['last_inconsistency_reason'][1]

            first_bad_vcf = metadata['vcfs'].index(bad_vcfs[0])
            if first_bad_vcf == nil
               raise "the metadata state is incoherent"
            end
            bad_samples = []
            metadata['samples'].each do |s| 
               if s['vcfid'] >= first_bad_vcf
                  bad_samples << s['name']
               end
            end

            number_of_good_samples = metadata['samples'].length - bad_samples.length
            number_of_good_vcfs = metadata['vcfs'].length - bad_vcfs.length

            update_operation = {
               '$push' => {
                  'IDs'     => {'$each' => [], '$slice'=> number_of_good_vcfs},
                  'QUALs'   => {'$each' => [], '$slice'=> number_of_good_vcfs}, 
                  'FILTERs' => {'$each' => [], '$slice'=> number_of_good_vcfs},
                  'INFOs'   => {'$each' => [], '$slice'=> number_of_good_vcfs},
                  'samples' => {'$each' => [], '$slice'=> number_of_good_samples}
                  }
               }

            db.collection(coll_name).update({}, update_operation, {:multi => true})

            metadata_update_operation = {
               '$set' => {'consistent' => true}, 
               '$push' => {
                  'vcfs'    => {'$each' => [], '$slice' => number_of_good_vcfs},
                  'headers' => {'$each' => [], '$slice' => number_of_good_vcfs},
                  'samples' => {'$each' => [], '$slice'=> number_of_good_samples}
                  },
               }
            db.collection('__METADATA__').update({'_id' => coll_name}, metadata_update_operation)
         else
            raise "unknown error state, fix_collection() doesn't know what to do"
         end
      end
   end
end