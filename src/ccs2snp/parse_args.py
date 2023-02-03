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
        "-v",
        "--version",
        action="version",
        version="%(prog)s {version}".format(version=program_version),
    )
    # subcommands: init
    subparsers = parser.add_subparsers(dest="sub", metavar="")

    # subcommands: call
    parser_call = subparsers.add_parser(
        "call",
        formatter_class=make_wide(argparse.ArgumentDefaultsHelpFormatter, w=180, h=60),
        help="detects somatic mutations from circular consensus seuqence (CCS) reads",
    )
    parser_call.add_argument(
        "-i",
        "--bam",
        type=str,
        required=True,
        help="minimap2 (parameters: -ax map-hifi --cs=short) aligned SAM/BAM files",
    )
    parser_call.add_argument(
        "--ref",
        type=str,
        required=True,
        help="reference FASTA file"
    )
    parser_call.add_argument(
        "--region",
        type=str,
        required=False,
        help="target chromosome",
    )
    parser_call.add_argument(
        "--region_list",
        type=str,
        required=False,
        help="list of target chromosomes separated by new line",
    )
    parser_call.add_argument(
        "--min_bq",
        type=int,
        default=20,
        required=False,
        help="minimum base quality score (BQ) threshold",
    )
    parser_call.add_argument(
        "--min_gq",
        type=int,
        default=20,
        required=False,
        help="minimum germline genotype quality (GQ) score ",
    )
    parser_call.add_argument(
        "--min_mapq",
        type=int,
        default=20,
        required=False,
        help="minimum mapping quality score",
    )
    parser_call.add_argument(
        "--min_ref_count",
        type=int,
        default=3,
        required=False,
        help="minimum reference allele depth at single base substitution site",
    )
    parser_call.add_argument(
        "--min_alt_count",
        type=int,
        default=1,
        required=False,
        help="minimum alternative allele depth at single base substitution site",
    )
    parser_call.add_argument(
        "--min_hap_count",
        type=int,
        default=3,
        required=False,
        help="minimum h0 and h1 haplotype count",
    )
    parser_call.add_argument(
        "--mismatch_window",
        type=int,
        default=20,
        required=False,
        help="mismatch window size",
    )
    parser_call.add_argument(
        "--max_mismatch_count",
        type=int,
        default=0,
        required=False,
        help="maximum number of mismatches within the mismatch window",
    )
    parser_call.add_argument(
        "--germline_snv_prior",
        type=float,
        default=1/(10**3),
        required=False,
        help="germline snv prior",
    )
    parser_call.add_argument(
        "--germline_indel_prior",
        type=float,
        default=1/(10**4),
        required=False,
        help="germline indel prior",
    )
    parser_call.add_argument(
        "-t",
        "--threads",
        type=int,
        default=1,
        required=False,
        help="number of threads to use",
    )
    parser_call.add_argument(
        "-o",
        "--out",
        type=str,
        required=True,
        help="VCF file to write the somatic substitutions",
    )
    # subcommands: sbs48
    parser_sbs48 = subparsers.add_parser(
        "sbs48",
        formatter_class=make_wide(argparse.ArgumentDefaultsHelpFormatter, w=180, h=60),
        help="returns sbs48 counts and barplot",
    )
    parser_sbs48.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help="himut VCF file to read somatic single base substitutions",
    )
    parser_sbs48.add_argument(
        "--ref", 
        type=str, 
        required=True, 
        help="reference FASTA file"
    )
    parser_sbs48.add_argument(
        "--region",
        type=str,
        required=False,
        help="target chromosome",
    )
    parser_sbs48.add_argument(
        "--region_list",
        type=str,
        required=False,
        help="list of target chromosomes separated by new line",
    )
    parser_sbs48.add_argument(
        "-o",
        "--output",
        type=str,
        required=True,
        help="file to return sbs48 counts (.tsv suffix)",
    )
    if len(arguments) == 0:
        parser.print_help()
        parser.exit()
    else:
        return parser, parser.parse_args(arguments)
