Gem::Specification.new do |s|
  s.name        = 'vcf-mongo'
  s.version     = '0.1.0'
  s.date        = '2014-08-14'
  s.summary     = "A small collection of scripts that uses MongoDB to store VCF files in collections."
  s.description = "A single collection can contain multiple VCF file that get soft-merged, meaning that you can conceptually think a collection as a big VCF file, but no information is lost, as some fields are 'multiplexed' (the `ID` field becomes an `IDs` array field, for example) leaving the choice of how to merge the data (when necessary) to the querying scripts (either provided or custom built by the user)."
  s.authors     = ["Loris Cro"]
  s.email       = 'l.cro@campus.unimib.it'
  s.files       = [
    "bin/vcf-import",
    "bin/vcf-admin",
    "bin/vcf-private",
    "lib/vcf-mongo.rb",
    "lib/vcf-mongo/admin.rb",
    "lib/vcf-mongo/import.rb",
    "htsjdk/sam-1.113.jar",
    "htsjdk/tribble-1.113.jar",
    "htsjdk/variant-1.113.jar"
  ]
  s.executables = [
    'vcf-import',
    'vcf-admin',
    'vcf-private'
  ]
  s.add_runtime_dependency "mongo"
  s.licenses       = ['BSD', 'MIT']
end