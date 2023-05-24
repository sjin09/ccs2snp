#!/usr/bin/env python3
__version__ = "0.0.1"
__author__ = "Sangjin Lee"

# modules
import sys
import himut.util
import ccs2snp.util
import ccs2snp.caller
from ccs2snp.parse_args import parse_args


def main():
    options = parse_args(program_version=__version__)
    himut.util.check_num_threads(options.threads)
    ccs2snp.caller.call_germline_snps(
        options.bam, # BAM file
        options.ref, # reference FASTA file
        options.region, # target contigs/scaffolds/chromosomes
        options.region_list, # target contigs/scaffolds/chromosomes fofn
        options.min_bq, # minimum base quality score: int
        options.min_gq, # minimum germline genotype quality score: int
        options.min_mapq, # int: 0-60
        options.min_ref_count, # number of reads supporting the reference base
        options.min_alt_count, # number of reads supporting the alterantive base
        options.mismatch_window, # mismatch window size
        options.max_mismatch_count, # maximum number of mismatches within a window
        options.germline_snv_prior,
        options.germline_indel_prior,
        options.threads,  # maxminum number of threads
        __version__,  # str
        options.out,  # output # ccs2snp vcf file
    )
    sys.exit(0)

if __name__ == "__main__":
    main()
