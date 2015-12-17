"""
For running locally and testing the conversion works.

@madprime: This takes ~4m on my laptop to convert a var file. Written to test
the package, but can also be used as a command-line conversion tool.
"""
import argparse
import os
import sys

from twobitreader import download as twobitdownload

import cgivar2vcf


def get_reference_genome_file(refseqdir, build):
    if not os.path.exists(refseqdir) or not os.path.isdir(refseqdir):
        raise ValueError("No directory at {}".format(refseqdir))
    twobitname = ''
    if build in ['b37', 'build 37', 'build37', '37', 'hg19']:
        twobitname = 'hg19.2bit'
    if not twobitname:
        raise ValueError('Genome bulid "{}" not supported.'.format(build))
    twobit_path = os.path.join(refseqdir, twobitname)
    if not os.path.exists(twobit_path):
        twobitdownload.save_genome('hg19', destdir=refseqdir)
    return twobit_path


def main(args):
    # Get local twobit file from its directory. Download and store if needed.
    twobit_path = get_reference_genome_file(args.refseqdir, build='b37')
    # Handle input
    if sys.stdin.isatty():  # false if data is piped in
        var_input = args.cgivarfile
    else:
        var_input = sys.stdin
    # Handle output
    if args.vcfoutfile:
        cgivar2vcf.convert_to_file(var_input,
                                   args.vcfoutfile,
                                   twobit_path)
    else:
        for line in cgivar2vcf.convert(var_input, twobit_path):
            print(line)


if __name__ == "__main__":
    # Parse options
    parser = argparse.ArgumentParser(
        description='Convert Complete Genomics var files to gVCF format.')
    parser.add_argument(
        '-d', '--refseqdir', metavar='REFSEQDIR', required=True,
        dest='refseqdir',
        help='Directory twobit reference genomes files are stored.')
    parser.add_argument(
        '-i', '--input', metavar='INPUTVARFILE',
        dest='cgivarfile',
        help='Path to Complete Genomics var file to convert. Can also be piped'
        ' in as standard input.')
    parser.add_argument(
        '-o', '--output', metavar='OUTPUTVCFFILE',
        dest='vcfoutfile',
        help='Path to where to save output VCF file.')
    args = parser.parse_args()
    main(args)
