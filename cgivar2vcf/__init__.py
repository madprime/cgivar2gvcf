#!/usr/bin/python
# Filename: cgivar2vcf.py
"""Conversion of Complete Genomics, Inc. (CGI) var files to VCF files."""
import bz2
import gzip
import re

import twobitreader


def auto_zip_open(filepath, mode):
    """Convenience function for opening potentially-compressed files."""
    if filepath.endswith('.gz'):
        outfile = gzip.open(filepath, mode)
    elif filepath.endswith('.bz2'):
        outfile = bz2.BZ2File(filepath, mode)
    else:
        outfile = open(filepath, mode)
    return outfile


def process_full_position(data, header):
    """Return genetic data when all alleles called on same line.

    Returns an array containing one item, a tuple of five items:
        (string) chromosome
        (string) start position (1-based)
        (array of strings) matching dbSNP entries
        (string) reference allele sequence
        (array of strings) the genome's allele sequences
    """
    feature_type = data[header['varType']]
    # Skip unmatchable, uncovered, or pseudoautosomal-in-X
    if (feature_type == 'no-ref' or feature_type.startswith('no-call') or
            feature_type.startswith('PAR-called-in-X')):
        return None
    # TODO: Don't skip REF, call using gVCF syntax.
    if feature_type == 'ref':
        return None

    # Skip low quality
    if ('varQuality' in header and
            data[header['varQuality']].startswith('VQLOW')):
        return None

    chrom = data[header['chromosome']]
    # One-based start to match VCF coordinates
    start = str(int(data[header['begin']]) + 1)
    ref_allele = data[header['reference']]
    alleles = [data[header['alleleSeq']]]
    dbsnp_data = []
    dbsnp_data = data[header['xRef']].split(';')
    return [(chrom, start, dbsnp_data, ref_allele, alleles)]


def process_allele(allele_data, dbsnp_data, header):
    """Combine data from multiple lines refering to a single allele.

    Returns three items in this order:
        (string) concatenated variant sequence (ie allele the genome has)
        (string) concatenated reference sequence
        (string) start position (1-based)
    """
    # One-based start to match VCF coordinates
    start = str(int(allele_data[0][header['begin']]) + 1)
    var_allele = ''
    ref_allele = ''
    for data in allele_data:
        # We reject allele data if any subset of the data has a no-call.
        if (data[header['varType']] == 'no-call' or
                ('varQuality' in header and
                 data[header['varQuality']].startswith('VQLOW'))):
            var_allele = None
            ref_allele = None
            break
        var_allele = var_allele + data[header['alleleSeq']]
        ref_allele = ref_allele + data[header['reference']]
        if data[header['xRef']]:
            for dbsnp_item in data[header['xRef']].split(';'):
                dbsnp_data.append(dbsnp_item.split(':')[1])
    return var_allele, ref_allele, start


def get_split_pos_lines(data, cgi_input, header):
    """Advance across split alleles and return data from each.

    CGI var file reports alleles separately for heterozygous sites:
    all variant or reference information is called for the first allele,
    then for the second. This function moves forward in the file to
    get lines for each (and ends up with one remaineder line as well).
    """
    s1_data = [data]
    s2_data = []
    next_data = cgi_input.next().rstrip('\n').split("\t")
    while next_data[header['allele']] == "1":
        s1_data.append(next_data)
        next_data = cgi_input.next().rstrip('\n').split("\t")
    while next_data[header['allele']] == "2":
        s2_data.append(next_data)
        next_data = cgi_input.next().rstrip('\n').split("\t")
    return s1_data, s2_data, next_data


def process_split_position(data, cgi_input, header):
    """Process CGI var where alleles are reported separately.

    Split positions report each allele with one or more lines. To ensure that
    we've read through all lines, we end up reading one line beyond.

    This function returns data for this position, then handles the remainder
    line by calling itself or process_full_position (as appropriate).

    Returns an array containing tuples with five items each:
        (string) chromosome
        (string) start position (1-based)
        (array of strings) matching dbSNP entries
        (string) reference allele sequence
        (array of strings) the genome's allele sequences
    """
    assert data[2] == "1"
    chrom = data[header['chromosome']]

    # Get all lines for each allele. Note that this means we'll end up with
    # data from one line ahead stored in 'next_data'; it will be handled at
    # the end.
    s1_data, s2_data, next_data = get_split_pos_lines(data, cgi_input, header)

    # Process all the lines to get concatenated sequences and other data.
    dbsnp_data = []
    a1_seq, ref_seq, start = process_allele(s1_data, dbsnp_data, header)
    a2_seq, r2_seq, a2_start = process_allele(s2_data, dbsnp_data, header)
    # clean dbsnp data
    dbsnp_data = [x for x in dbsnp_data if x]
    if (a1_seq or ref_seq) and (a2_seq or r2_seq):
        # Check that reference sequence and positions match.
        assert ref_seq == r2_seq
        assert start == a2_start
        yield (chrom, start, dbsnp_data, ref_seq, [a1_seq, a2_seq])

    # Handle the remaining line (may recursively call this function if it's
    # the start of a new region with separated allele calls).
    if next_data[2] == "all" or next_data[1] == "1":
        out = process_full_position(next_data, header)
    else:
        out = process_split_position(next_data, cgi_input, header)
    if out:
        for entry in out:
            yield entry


def vcf_line(input_data, reference):
    """
    Convert the var files information into VCF format.

    This is nontrivial because the var file can contain zero-length variants,
    which is not allowed by VCF. To handle these cases, we "move backwards"
    by one position, look up the reference sequence, and add that.

    The returned line is a very simple, VCF-valid row containing the
    genome's data for this position.
    """
    chrom, start, dbsnp_data, ref_allele, genome_alleles = input_data

    # VCF start position is one less than var file's indexing.
    start = int(start) - 1

    # VCF doesn't allow zero-length sequences. If we have this situation,
    # move the start backwards by one position, get that reference base,
    # and prepend this base to all sequences.
    if len(ref_allele) == 0 or 0 in [len(v) for v in genome_alleles]:
        start = start - 1
        prepend = reference[chrom][start]
        ref_allele = prepend + ref_allele
        genome_alleles = [prepend + v for v in genome_alleles]

    # Figure out what our alternate alleles are.
    alt_alleles = []
    for allele in genome_alleles:
        if allele not in [ref_allele] + alt_alleles:
            alt_alleles.append(allele)

    # Combine ref and alt for the full set of alleles, used for indexing.
    alleles = [ref_allele] + alt_alleles

    # Get the indexed genotype.
    genotype = '/'.join([str(alleles.index(x)) for x in genome_alleles])

    # Get dbSNP IDs.
    dbsnp_cleaned = []
    for dbsnp in dbsnp_data:
        if dbsnp not in dbsnp_cleaned:
            dbsnp_cleaned.append(dbsnp)
    if dbsnp_cleaned:
        id_field = ';'.join(dbsnp_cleaned)
    else:
        id_field = '.'

    return '\t'.join([chrom, str(start), id_field, ref_allele,
                      ','.join(alt_alleles), '.', '.', '.',
                      'GT', genotype])


def process_next_position(data, cgi_data, header, reference):
    """
    Determine appropriate processing to get data, then convert it to VCF

    There are two types of lines in the var file:
    - "full position": single allele (hemizygous) or all-allele line
        All alleles at this position are represented in this line.
        This is handled with "process_full_position".
    - "split position": each of two alleles is reported separately. There will
        be at least two lines, one for each allele (but potentially more).
        This is handled with "process_split_position".

    Because the number of lines used for separately reported alleles is
    unknown, process_split_position will always read ahead to the next
    "full position" and return that as well.

    So the returned line formats are consistent, process_next_position
    returns an array, even if there's only one line.
    """
    if data[2] == "all" or data[1] == "1":
        # The output from process_full_position is an array, so it can be
        # treated in the same manner as process_split_position output.
        out = process_full_position(data, header)
    else:
        assert data[2] == "1"
        # The output from process_split_position is a generator, and may end
        # up calling itself recursively.
        out = process_split_position(data, cgi_data, header)
    if out:
        return [vcf_line(l, reference) for l in out]


def convert(cgi_data, twobit_ref):
    """Generator that converts CGI var data to VCF-formated strings"""

    # Set up CGI input. Default is to assume a str generator.
    if isinstance(cgi_data, basestring):
        cgi_data = auto_zip_open(cgi_data, 'rb')

    # Set up TwoBitFile for retrieving reference sequences.
    reference = twobitreader.TwoBitFile(twobit_ref)

    for line in cgi_data:
        # Skip header lines.
        if re.search(r'^\W*$', line) or line.startswith('#'):
            continue

        # Store header row labels.
        if line.startswith('>'):
            header_data = line.lstrip('>').rstrip('\n').split('\t')
            header = {header_data[i]: i for i in range(len(header_data))}
            continue

        # If we reach this point, this is a line that contains data.
        data = line.rstrip('\n').split("\t")

        out = process_next_position(data, cgi_data, header, reference)

        # process_next_position returns an array of one or more lines
        if out:
            for line in out:
                yield line


def convert_to_file(cgi_input, output_file, twobit_ref):
    """Convert a CGI var file and output VCF-formatted data to file"""

    if isinstance(output_file, basestring):
        output_file = auto_zip_open(output_file, 'wb')

    conversion = convert(cgi_input, twobit_ref)  # set up generator
    for line in conversion:
        output_file.write(line + "\n")
    output_file.close()
