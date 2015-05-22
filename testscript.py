"""
For testing the conversion works. Takes ~6m on my laptop to run. -mpball
"""

import sys

from optparse import OptionParser

import cgivar2vcf


def main():
    # Parse options
    usage = "\n%prog -i inputfile [-o outputfile]\n" \
            + "%prog [-o outputfile] < inputfile"
    parser = OptionParser(usage=usage)
    parser.add_option("-i", "--input", dest="inputfile",
                      help="read CGI data from INFILE (uncompressed, .gz," +
                      " or .bz2)", metavar="INFILE")
    parser.add_option("-o", "--output", dest="outputfile",
                      help="write report to OUTFILE (uncompressed, " +
                      "*.gz, or *.bz2)", metavar="OUTFILE")
    parser.add_option("-r", "--ref2bit", dest="twobitref",
                      help="2bit reference genome file",
                      metavar="TWOBITREF")
    options, _ = parser.parse_args()
    # Handle input
    if sys.stdin.isatty():  # false if data is piped in
        var_input = options.inputfile
    else:
        var_input = sys.stdin
    # Handle output
    if options.outputfile:
        cgivar2vcf.convert_to_file(var_input,
                                   options.outputfile,
                                   options.twobitref)
    else:
        for line in cgivar2vcf.convert(var_input, options.twobitref):
            print line

if __name__ == "__main__":
    main()
