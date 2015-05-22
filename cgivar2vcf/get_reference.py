from __future__ import absolute_import
from .twobit import TwoBitFile


def get_reference_allele(chrom, start, hg19_path):
    twobit_file = TwoBitFile(hg19_path)
    end = start + 1
    try:
        refallele = twobit_file[chrom][start:end]
    except TypeError:
        refallele = 'N'
    return refallele
