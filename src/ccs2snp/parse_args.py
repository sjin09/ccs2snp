# modules
import sys
import warnings
import argparse


def make_wide(formatter, w=120, h=36):
    """Return a wider HelpFormatter, if possible."""
    try:
        # https://stackoverflow.com/a/5464440
        # beware: "Only the name of this class is considered a public API."
        kwargs = {"width": w, "max_help_position": h}
        formatter(None, **kwargs)
        return lambda prog: formatter(prog, **kwargs)
    except TypeError:
        warnings.warn("argparse help formatter failed, falling back.")
        return formatter


def parse_args(program_version, arguments=sys.argv[1:]):
    # main_arguments
    parser = argparse.ArgumentParser(
        add_help=True,
        formatter_class=make_wide(argparse.ArgumentDefaultsHelpFormatter),
        description="ccs2snp identifies germline SNPs from PacBio CCS reads",
    )
    parser.add_argument(
        "-i",
        "--bam",
        type=str,
        required=True,
        help="minimap2 (parameters: -ax map-hifi --cs=short) aligned SAM/BAM files",
    )
    parser.add_argument(
        "--ref",
        type=str,
        required=True,
        help="reference FASTA file"
    )
    parser.add_argument(
        "--region",
        type=str,
        required=False,
        help="target chromosome",
    )
    parser.add_argument(
        "--region_list",
        type=str,
        required=False,
        help="list of target chromosomes separated by new line",
    )
    parser.add_argument(
        "--min_bq",
        type=int,
        default=20,
        required=False,
        help="minimum base quality score (BQ) threshold",
    )
    parser.add_argument(
        "--min_gq",
        type=int,
        default=20,
        required=False,
        help="minimum germline genotype quality (GQ) score ",
    )
    parser.add_argument(
        "--min_mapq",
        type=int,
        default=20,
        required=False,
        help="minimum mapping quality score",
    )
    parser.add_argument(
        "--min_ref_count",
        type=int,
        default=3,
        required=False,
        help="minimum reference allele depth at single base substitution site",
    )
    parser.add_argument(
        "--min_alt_count",
        type=int,
        default=1,
        required=False,
        help="minimum alternative allele depth at single base substitution site",
    )
    parser.add_argument(
        "--mismatch_window",
        type=int,
        default=20,
        required=False,
        help="mismatch window size",
    )
    parser.add_argument(
        "--max_mismatch_count",
        type=int,
        default=0,
        required=False,
        help="maximum number of mismatches within the mismatch window",
    )
    parser.add_argument(
        "--germline_snv_prior",
        type=float,
        default=1/(10**3),
        required=False,
        help="germline snv prior",
    )
    parser.add_argument(
        "--germline_indel_prior",
        type=float,
        default=1/(10**4),
        required=False,
        help="germline indel prior",
    )
    parser.add_argument(
        "-t",
        "--threads",
        type=int,
        default=1,
        required=False,
        help="number of threads to use",
    )
    parser.add_argument(
        "-o",
        "--out",
        type=str,
        required=True,
        help="VCF file to write the somatic substitutions",
    )
    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version="%(prog)s {version}".format(version=program_version),
    )
    if len(arguments) == 0:
        parser.print_help()
        parser.exit()
    else:
        return parser.parse_args(arguments)
