#!/usr/bin/env python3
__version__ = "1.0.0"
__author__ = "Sangjin Lee"

# modules

import himut.util
import ccs2snp.util
import ccs2snp.caller
import ccs2snp.mutlib
from ccs2snp.parse_args import parse_args


def main():
    parser, options = parse_args(program_version=__version__)
    if options.sub == "call":  # call single nucleotide polymorphisms
        himut.util.check_num_threads(options.threads)
        ccs2snp.caller.call_germline_snps(
            options.bam, # BAM file
            options.ref, # reference FASTA file
            options.region, # target contigs/scaffolds/chromosomes
            options.region_list, # target contigs/scaffolds/chromosomes fofn
            options.min_bq, # minimum basequality score: int
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
    elif options.sub == "sbs48":  # returns sbs96 counts
        sample = ccs2snp.vcflib.get_sample(options.input)
        ccs2snp.util.check_mutpatterns_input_exists(
            options.input,
            options.ref,
            options.region,
            options.region_list,
            options.output,
        )
        ccs2snp.mutlib.dump_sbs48_counts(
            options.input,
            options.ref,
            options.region,
            options.region_list,
            options.output,
        )
        ccs2snp.mutlib.dump_sbs48_plt(
            options.output, sample, "{}.pdf".format(options.output)
        )
    else:
        print("The subcommand does not exist!\n")
        parser.print_help()
        parser.exit()


if __name__ == "__main__":
    main()
    ccs2snp.util.exit()
