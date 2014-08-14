vcf-mongo
=========

A small collection of scripts that uses MongoDB to store VCF files in collections.

A single collection can contain multiple VCF file that get soft-merged, meaning that you can conceptually think a collection as a big VCF file, but no information is lost, as some fields are 'multiplexed' (the `ID` field becomes an `IDs` array field, for example) leaving the choice of how to merge the data (when necessary) to the querying scripts (either provided or custom built by the user).

A soft-merge approach should offer a somewhat natural way of querying the database while still offering good flexibility in terms of fruition.

An example script that queries for private variants is also provided.

REQUIREMENTS
------------

* A running MongoDB instance ([http://mongodb.org](http://mongodb.org))
* A working JRuby setup([http://jruby.org/](http://jruby.org/))
* The MongoDB driver for JRuby (`jgem install mongo`)

USAGE
-----

### Importing VCF files ###

Basic example:

`jruby vcf-import.rb myCollection file1.vcf file2.vcf.gz ...`

You can use the `-h` flag to see all options.

Parametric options:
* Connection: `--address`, `--port`
* Database selection: `--database`
* Number of threads that do the merge operations / execute the import queries: `--mongo-threads`, `--merge-threads`
* Size of each thread's work buffer: `--parser-buffer-size`, `--merger-buffer-size`, `--mongo-buffer-size`

Switches:
* Enable the (slower) import into an existing collection: `--append`
* Drop bad records without causing the whole opration to fail (bad as in 'rejected by the HTSJDK parser'): `--drop-bad-records` (a warning is printed to STDERR for each dropped record)
* Disable the showing of import speed and buffer saturation (which would take a lot of space when logging STDOUT): `--no-progress`

Both in case of a new database or a new collection, the script does the initialization for you, you don't need to use MongoDB directly for any of these operations.

The import operations are done via the following pipeline:

* Each VCF file summons a parser in a different thread.
* Each parser-thread yields Record objects to the Aligner thread.
* Multiple merge-threads mix the corresponding records (same CHROM:POS) yielded from the Aligner thread.
* Multiple mongo-threads perform the import of these multi-records via insert (or upsert, in case of an `--append` operation).

The best combination of buffer size and thread number depends on your machine, wether the MongoDB instance is local or remote, how many VCF files you're importing and how complex each file is (multisample, lots of subfields).

The script also provides a count of how many items are inside each work queue. You can use that information to get an insight on where the bottleneck is.

As of now there are two *big-ish* things that need to be fixed:
* There is no completion percentage of the operation because the HTSJDK parser hides the InputStream instance (so we can't `tell()` the current position): [issue #63](https://github.com/samtools/htsjdk/issues/63).
* The actual *hard* work of parsing all subfields is currently done by the parser threads. By design these kind of operations should be perfomed by the merger threads but apparently multiple Record instances yielded by the same parser (these records are called `VariantContext` by the HTSJDK parser) are not thread-safe (somewhere there is a pesky shared resource). This makes the import somewhat slower (especially in cases of a single multisample file) than is could be because records from the same VCF file have to be parsed in a serial fashon.

There are other minor things that could be noted but these two have an actual impact on the script's usage. For additional information (more from a development point of view than an usage one), check out the wiki section of this repo.


### Managing your collections ###

Use `vcf-admin.rb`. 

The commands are: 
* `list`: to obtain the list of all collection present or details about a single collection
* `rename`, `delete`: what the command says
* `check`: collections that have an import operation that didn't complete are flagged as *inconsistent*; this command lists them
* `fix`: fix an inconsistent collection

You can use the `-h` flag to see all options.

Destructive operations (`delete`, `fix`) require confirmation unless the `--force` option is specified.


##### Checking and fixing collections ####
If an import operation fails the target collection is marked as *inconsistent*.
If the import operation was into an already existing collection (`--append`), the operation can be reverted, otherwise fixing the collection means deleting it (as the possibility of resuming an import is not yet available).

### Querying your collections ###

TODO



CONTRIBUTING
------------

TODO







