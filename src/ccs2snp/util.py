import os
import sys
import pysam
import pyfastx
from typing import Dict, List, Tuple


def is_bam_file_corrupt(bam_file: str, chrom_lst: List[str]):

    hsh = {}
    alignments = pysam.AlignmentFile(bam_file, "rb")
    for chrom in chrom_lst:
        hsh[chrom] = alignments.count(chrom)

    state = 0
    for chrom in chrom_lst:
        if hsh[chrom] == 0:
            print(
                "bam_file: {} does not have read alignments for chrom: {}".format(
                    bam_file, chrom
                )
            )
            state = 1

    if state == 1:
        print("--bam {} might be corrupted".format(bam_file))
        return 1
    else:
        return 0


def is_ref_file_corrupt(ref_file: str, chrom_lst: List[str]):
    state = 0
    refseq = pyfastx.Fasta(ref_file)
    for chrom in chrom_lst:
        if chrom not in refseq:
            print("chrom :{} is missing from --ref {}".format(chrom, ref_file))
            state = 1

    if state == 1:
        print("--ref {} might be corrupted".format(ref_file))
        return 1
    else:
        return 0


def is_ref_file(ref_file: str):

    if ref_file is None:
        print("Please provide a path to the reference FASTA file")
        return 1
    else:
        if os.path.exists(ref_file):
            if ref_file.endswith((".fa", ".fa.gz", ".fasta", ".fasta.gz")):
                if os.path.getsize(ref_file) == 0:
                    return 1
                else:
                    return 0
            else:
                print("Did you provide a FASTA file?")
                return 1
        else:
            print("reference FASTA file is missing")
            return 1



def is_bam_file(bam_file: str):

    if bam_file is None:
        print("Please provide the path to the BAM file")
        return 1
    else:
        if bam_file.endswith(".bam"):
            idxfile = "{}.bai".format(bam_file)
            if os.path.exists(bam_file) and os.path.exists(idxfile):
                if os.path.getsize(bam_file) != 0 and os.path.getsize(idxfile) != 0:
                    return 0
                else:
                    return 1
            elif not os.path.exists(bam_file) and os.path.exists(idxfile):
                print("BAM file is missing")
                return 1
            elif os.path.exists(bam_file) and not os.path.exists(idxfile):
                print("BAM index file is missing")
                print("Use samtools index to index your BAM file")
                return 1
            elif not os.path.exists(bam_file) and not os.path.exists(idxfile):
                print("BAM file is missing")
                print("BAM index file is missing")
                return 1
        else:
            print("Did you provide a BAM file?")
            print("BAM files have to a .bam suffix")
            return 1


def is_out_file(out_file: str):
    if out_file is None:
        return 1
    else:
        if out_file.endswith(".vcf"):
            return 0
        else:
            print("ccs2snp doesn't recognise the suffix of the file")
            return 1


def is_input_exist(
    bam_file: str,
    ref_file: str,
    out_file: str,
):

    counter = 0
    counter += is_ref_file(ref_file)
    counter += is_bam_file(bam_file)
    counter += is_out_file(out_file)
    if counter > 0:
        print("One or more inputs and parameters are missing")
        print("Please provide the correct inputs and parameters")
        sys.exit(0)


def is_input_corrupt(
    bam_file: str,
    ref_file: str,
    chrom_lst: List[str]
):

    counter = 0
    counter += is_ref_file_corrupt(ref_file, chrom_lst)
    counter += is_bam_file_corrupt(bam_file, chrom_lst)

