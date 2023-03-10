import os
import sys
import tabix
import pysam
import cyvcf2
import psutil
import pyfastx
import natsort
import himut.vcflib
from collections import defaultdict
from typing import List, Dict, Tuple




def exit():
    print("exiting ccs2snp")
    sys.exit(0)





def check_bam_file(
    bam_file: str,
    chrom_lst: List[str],
):

    if bam_file is None:
        print("Please provide the path to the BAM file")
        return 1
    else:
        if bam_file.endswith(".bam"):
            idxfile = "{}.bai".format(bam_file)
            if os.path.exists(bam_file) and os.path.exists(idxfile):
                if os.path.getsize(bam_file) != 0 and os.path.getsize(idxfile) != 0:
                    if is_bam_file_corrupt(bam_file, chrom_lst):
                        return 0
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


def is_vcf_file_corrupt(
    param: str, vcf_file: str, chrom_lst: List[str], tname2tsize: Dict[str, int]
) -> bool:

    hsh = defaultdict(lambda: 0)
    if vcf_file.endswith(".vcf"):
        for line in open(vcf_file):
            if line.startswith("#"):
                continue
            arr = line.strip().split()
            chrom = arr[0]
            ref = arr[3]
            alt_lst = arr[4].split(",")
            if arr[6] == "PASS" and len(alt_lst) == 1:
                alt = alt_lst[0]
                if len(ref) == 1 and len(alt) == 1:
                    hsh[chrom] += 1
    elif vcf_file.endswith(".bgz"):
        tb = tabix.open(vcf_file)
        for chrom in chrom_lst:
            try:
                records = tb.query(chrom, 0, tname2tsize[chrom])
                hsh[chrom] = len(list(records))
            except tabix.TabixError:
                continue

    state = 0
    for chrom in chrom_lst:
        if hsh[chrom] == 0:
            print(
                "{} {} does not have any variants on chrom: {}".format(
                    param, vcf_file, chrom
                )
            )
            state = 1
    if state == 1:
        print("{} {} might be corrupted".format(param, vcf_file))
        return 1
    else:
        return 0


def check_vcf_file(
    param: str, vcf_file: str, chrom_lst: List[str], tname2tsize: Dict[str, int]
):
    if vcf_file is None:
        print(
            "Please provide the path to the VCF file for the following argument: {}".format(
                param
            )
        )
        return 1
    else:
        if os.path.exists(vcf_file):
            if vcf_file.endswith(".vcf"):
                if os.path.getsize(vcf_file) == 0:
                    return 1
                else:
                    if is_vcf_file_corrupt(param, vcf_file, chrom_lst, tname2tsize):
                        return 0
                    return 0
            elif vcf_file.endswith(".bgz"):
                tbi_file = vcf_file + ".tbi"
                if os.path.exists(tbi_file):
                    if is_vcf_file_corrupt(param, vcf_file, chrom_lst, tname2tsize):
                        return 0
                    return 0
                else:
                    print(
                        "tabix index file does not exist {} {}".format(param, vcf_file)
                    )
                    return 1
            elif vcf_file.endswith(".gz"):
                print(
                    "himut doesn't support loading of gzip compressed VCF files {} {}".format(
                        param, vcf_file
                    )
                )
                return 1
            else:
                print(
                    "VCF file must have the .vcf suffix {} {}".format(param, vcf_file)
                )
                return 1
        else:
            print("{} {} file is missing".format(param, vcf_file))
            return 1


def check_out_file(out_file: str):
    if out_file is None:
        return 1
    else:
        if out_file.endswith(".vcf"):
            return 0
        elif out_file.endswith(".gz"):
            print("himut doesn't support return gzip compressed VCF files")
            return 1
        elif out_file.endswith(".bgz"):
            print("himut doesn't support return bgzip compressed VCF files")
            return 1
        else:
            print("himut doesn't recognise the suffix of the file")
            return 1


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


def check_ref_file(ref_file: str, chrom_lst: List[str]) -> int:

    if ref_file is None:
        print("Please provide the path to the reference FASTA file")
        return 1
    else:
        if os.path.exists(ref_file):
            if ref_file.endswith((".fa", ".fa.gz", ".fasta", ".fasta.gz")):
                if os.path.getsize(ref_file) == 0:
                    return 1
                else:
                    if is_ref_file_corrupt(ref_file, chrom_lst):
                        return 1
                    else:
                        return 0
            else:
                print("Did you provide a FASTA file?")
                return 1
        else:
            print("reference FASTA file is missing")
            return 1


def check_phased_vcf_file(
    vcf_file: str, chrom_lst: List[str], tname2tsize: Dict[str, int]
):
    if vcf_file is None:
        print(
            "Please provide the path to the VCF file for the following argument: --phased_vcf"
        )
        return 1
    else:
        if os.path.exists(vcf_file):
            if vcf_file.endswith(".vcf"):
                if os.path.getsize(vcf_file) == 0:
                    return 1
                else:
                    vcf_header_lst = []
                    for line in open(vcf_file).readlines():
                        if line.startswith("#"):
                            vcf_header_lst.append(line.strip())
                        if line.startswith("#CHROM"):
                            break

                    if any(
                        "##FORMAT=<ID=PS" in line.strip() for line in vcf_header_lst
                    ):
                        if is_vcf_file_corrupt(
                            "--phased_vcf", vcf_file, chrom_lst, tname2tsize
                        ):
                            return 1
                        else:
                            return 0
                    else:
                        print("VCF file is not phased")
                        return 1

            elif vcf_file.endswith(".bgz"):
                tbi_file = vcf_file + ".tbi"
                if os.path.exists(tbi_file):
                    v = cyvcf2.VCF(vcf_file)
                    vcf_header_lst = v.raw_header.strip().split()
                    if any(
                        "##FORMAT=<ID=PS" in line.strip() for line in vcf_header_lst
                    ):
                        if is_vcf_file_corrupt(
                            "--phased_vcf", vcf_file, chrom_lst, tname2tsize
                        ):
                            return 1
                        else:
                            return 0
                    else:
                        print("VCF file is not phased")
                        return 1
                else:
                    print("tabix index file does not exist")
                    return 1
            elif vcf_file.endswith(".gz"):
                print("himut doesn't support loading of gzip compressed VCF files")
                return 1
            else:
                print("VCF file must have the .vcf suffix")
                return 1
        else:
            print("{} file is missing".format(vcf_file))
            return 1


def check_caller_input_exists(
    bam_file: str,
    ref_file: str,
    vcf_file: str,
    phased_vcf_file: str,
    common_snps: str,
    panel_of_normals: str,
    chrom_lst: List[str],
    tname2tsize: Dict[str, int],
    phase: bool,
    non_human_sample: bool,
    create_panel_of_normals: bool,
    out_file: str,
) -> None:

    counter = 0
    counter += check_bam_file(bam_file, chrom_lst)
    counter += check_out_file(out_file)
    if phase:
        counter += check_phased_vcf_file(phased_vcf_file, chrom_lst, tname2tsize)

    if non_human_sample:
        counter += check_ref_file(ref_file, chrom_lst)
        counter += check_vcf_file("--vcf", vcf_file, chrom_lst, tname2tsize)

    if not non_human_sample and create_panel_of_normals:
        counter += check_vcf_file("--common_snps", common_snps, chrom_lst, tname2tsize)

    elif not non_human_sample and not create_panel_of_normals:
        counter += check_vcf_file("--common_snps", common_snps, chrom_lst, tname2tsize)
        counter += check_vcf_file(
            "--panel_of_normals", panel_of_normals, chrom_lst, tname2tsize
        )

    if counter > 0:
        print("One or more inputs and parameters are missing")
        print("Please provide the correct inputs and parameters")
        exit()


def check_phaser_input_exists(
    bam_file: str,
    vcf_file: str,
    out_file: str,
    chrom_lst: List[str],
    tname2tsize: Dict[str, int],
) -> None:

    counter = 0
    counter += check_bam_file(bam_file, chrom_lst)
    counter += check_vcf_file("--vcf", vcf_file, chrom_lst, tname2tsize)
    counter += check_out_file(out_file)
    if not out_file.endswith(".phased.vcf"):
        print("Please use the suffix .phased.vcf for the output file")
        counter += 1

    if counter > 0:
        print("One or more inputs and parameters are missing")
        print("Please provide the correct inputs and parameters")
        exit()


def check_mutpatterns_input_exists(
    vcf_file: str,
    ref_file: str,
    region: str,
    region_list: str,
    tname2tsize: Dict[str, int],
    out_file: str,
):

    counter = 0
    if region is None and region_list is not None:
        print(
            "himut will return SBS96/DBS78 counts from chromosomes and contings in {} file".format(
                region_list
            )
        )
    elif region is not None and region_list is None:
        print("himut will return SBS96/DBS78 counts from {}".format(region))
    elif region is not None and region_list is not None:
        print(
            "Please provide input for --region or --region_list parameter and not for both parameters"
        )
        counter = 1
    elif region is None and region_list is None:
        print("himut will return SBS96/DBS78 counts from all chromosomes and contigs")

    chrom_lst, _ = load_loci(region, region_list, tname2tsize)
    counter += check_vcf_file("--input", vcf_file, chrom_lst)
    counter += check_ref_file(ref_file, chrom_lst)
    if out_file is None:
        counter += 1
    if counter > 0:
        print("One or more inputs and parameters are missing")
        print("Please provide the correct inputs and parameters")
        exit()


def check_normcounts_input_exists(
    bam_file: str,
    ref_file: str,
    sbs_file: str,
    vcf_file: str,
    phased_vcf_file: str,
    common_snps: str,
    panel_of_normals: str,
    chrom_lst: List[str],
    tname2tsize: Dict[str, int],
    phase: bool,
    reference_sample: bool,
    non_human_sample: bool,
    out_file: str,
) -> None:

    counter = 0
    counter += check_bam_file(bam_file, chrom_lst) 
    counter += check_ref_file(ref_file, chrom_lst)
    counter += check_vcf_file("--sbs", sbs_file, chrom_lst, tname2tsize)
    
    if phase:
        counter += check_phased_vcf_file(phased_vcf_file, chrom_lst, tname2tsize)

    if reference_sample:
        counter += check_vcf_file("--vcf", vcf_file, chrom_lst, tname2tsize)

    if not non_human_sample:
        counter += check_vcf_file("--common_snps", common_snps, chrom_lst, tname2tsize)
        counter += check_vcf_file(
            "--panel_of_normals", panel_of_normals, chrom_lst, tname2tsize
        )

    if out_file is None:
        counter += 1
    if counter > 0:
        print("One or more inputs and parameters are missing")
        print("Please provide the correct inputs and parameters")
        exit()


def get_truncated_float(f: float) -> float:

    flst = [round(f, i) for i in range(1, 10)]
    mindex = max([j for j, k in enumerate(flst) if k == 0.0]) + 1
    tf = flst[mindex]
    return tf
