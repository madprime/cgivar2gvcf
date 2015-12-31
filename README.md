# cgivar2gvcf
Conversion of Complete Genomics var file to gVCF

This conversion assumes that the genome file is build 37 and does not currently
support other builds.

## Installation

This is Python package in the Python Package Index (PyPI). You can install it with
`pip install cgivar2gvcf`. This installation will also install the
`twobitreader` package.

## Command line usage

You can run this tool on the command line like this:

`python -m cgivar2gvcf -h`

The above command will display the program's options.

Notably, you need a copy of the UCSC 2bit reference genome to perform conversion.
The command line tool expects you to provide a directory where this file exists
(it should have the name `hg19.2bit`). If it's not present, the tool with download
a copy into this directory.

An example command for a variant-only VCF file (not gVCF):

`python -m cgivar2gvcf -d files/ -i var-GS00253-ASM.tsv.bz2 --var-only -o GS00253-vcf-from-var.vcf.bz2`

## Module usage

### convert_to_file

Writes the new VCF file to the specified output destination.

`convert_to_file(cgi_input, output_file, twobit_ref, twobit_name, var_only=False)`

* **cgi_input** - Inputed Complete Genomics var file. Can be a file-like object (i.e. already open) or a path to a file (uncompressed, gzip, or bzip2 compressed).
* **output_file** - Where to put the result. A "gz" suffix results in gzip compression, a "bz2" suffix results in bzip2 compression.
* **twobit_ref** - Path to the UCSC twobit reference genome file
* **twobit_name** - name of the twobit reference file (so we can cite this in the VCF header as the "reference")
* **var_only** - (optional) default false. Set to true if you only want variant lines.

### convert

Returns a generator object that yields lines of the VCF file.

`convert(cgi_input, twobit_ref, twobit_name, var_only=False)`

* **cgi_input** - Inputed Complete Genomics var file. Can be a file-like object (i.e. already open) or a path to a file (uncompressed, gzip, or bzip2 compressed).
* **twobit_ref** - Path to the UCSC twobit reference genome file
* **twobit_name** - name of the twobit reference file (so we can cite this in the VCF header as the "reference")
* **var_only** - (optional) default false. Set to true if you only want variant lines.

### get_reference_genome_file

Convenience function for finding the UCSC reference genome in a specified directory,
and downloading it if it's not present.

`get_reference_genome_file(refseqdir, build)`

* **refseqdir** - The directory that should have the 2bit file (it will be downloaded to this directory if not already present)
* **build** - The genome build. There is currently only one option: build 37. (You can use 'b37'.) It's here so we can be explicit about this.
