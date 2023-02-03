import tabix
import cyvcf2
import pyfastx
import natsort
import numpy as np
import himut.util
import himut.bamlib
from datetime import datetime
from collections import defaultdict
from typing import Dict, List, Set, Tuple


def get_vcf_header(
    bam_file: str,
    ref_file: str,
    region: str,
    region_list: str,
    tname2tsize: Dict[str, int],
    min_bq: int,
    min_gq: int,
    min_mapq: int,
    md_threshold: int,
    min_ref_count: int,
    min_alt_count: int,
    mismatch_window: int,
    max_mismatch_count: int,
    germline_snv_prior: float,
    germline_indel_prior: float,
    threads: int,
    version: str,
    out_file: str,
) -> str:

    vcf_header_lst = [
        "##fileformat=VCFv4.2",
        "##fileDate={}".format(datetime.now().strftime("%d%m%Y")),
        "##source=ccs2snp",
        "##source_version={}".format(version),
        '##FILTER=<ID=PASS,Description="All filters passed">',
        '##FILTER=<ID=LowBQ,Description="Base quality score is below the minimum base quality score of {}">'.format(
            min_bq
        ),
        '##FILTER=<ID=LowGQ,Description="Germline genotype quality score is below the minimum genotype quality score of {}">'.format(
            min_gq
        ),
        '##FILTER=<ID=IndelSite,Description="Somatic substitution at indel site is not considered">',
        '##FILTER=<ID=HetSite,Description="Somatic substitution at heterzygous SNP site is not considered">',
        '##FILTER=<ID=HetAltSite,Description="Somatic substitution at tri-allelic SNP site is not considered">',
        '##FILTER=<ID=HomAltSite,Description="Somatic substitution at homozygous alternative SNP site is not considered">',
        '##FILTER=<ID=LowDepth,Description="Read depth is below the minimum reference allele and/or alterantive allele depth threshold">',
        '##FILTER=<ID=HighDepth,Description="Read depth is above the maximum depth threshold of {:.1f}">'.format(
            md_threshold
        ),
        '##FILTER=<ID=MismatchConflict,Description="Substitution is found next to a mismatch within a given mismatch window">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '##FORMAT=<ID=GQ,Number=1,Type=String,Description="Genotype quality score">',
        '##FORMAT=<ID=BQ,Number=1,Type=Float,Description="Base quality score">',
        '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">',
        '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Read depth for each allele">',
        '##FORMAT=<ID=VAF,Number=A,Type=Float,Description="Variant allele fractions">',
        '##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set">',
    ]
    for tname in natsort.natsorted(list(tname2tsize.keys())):
        vcf_header_lst.append(
            "##contig=<ID={},length={}>".format(tname, tname2tsize[tname])
        )

    if region is None and region_list is not None:
        region_param = "--region_list {}".format(region_list)
    elif region is not None and region_list is None:
        region_param = "--region {}".format(region)
    elif region is not None and region_list is not None:
        region_param = "--region_list {}".format(region_list)

    cmdline = "##himut_command=himut call -i {} --ref {} {} --min_bq {} --min_gq {} --min_mapq {} --min_ref_count {} --min_alt_count {} --mismatch_window {} --max_mismatch_count {} --germline_snv_prior {} --germline_indel_prior {} --threads {} -o {}".format(
        bam_file,
        ref_file,
        region_param,
        min_bq,
        min_gq,
        min_mapq,
        min_ref_count,
        min_alt_count,
        mismatch_window,
        max_mismatch_count,
        germline_snv_prior,
        germline_indel_prior,
        threads,
        out_file,
    )
    vcf_header_lst.append(cmdline)
    vcf_header_lst.append(
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}".format(
            himut.bamlib.get_sample(bam_file)
        )
    )
    vcf_header = "\n".join(vcf_header_lst)
    return vcf_header



def dump_snp(
    vcf_file: str,
    vcf_header: str,
    chrom_lst: List[str],
    chrom2snp_lst: Dict[
        str, List[Tuple[str, int, str, str, str, str, int, int, int, float]]
    ],
):

    if not vcf_file.endswith(".vcf"):
        print("VCF file must have .vcf suffix")
        himut.util.exit()

    o = open(vcf_file, "w")
    o.write("{}\n".format(vcf_header))
    for chrom in chrom_lst:
        for (
            chrom,
            pos,
            ref,
            alt,
            state,
            vgt,
            gq,
            bq,
            read_depth,
            ref_count,
            alt_count,
            vaf,
        ) in chrom2snp_lst[chrom]:
            if vgt == "1/2":
                o.write(
                    "{}\t{}\t.\t{}\t{}\t.\t{}\t.\tGT:GQ:BQ:DP:AD:VAF\t{}:{}:{}:{:0.0f}:{:0.0f},{}:{}\n".format(
                        chrom,
                        pos,
                        ref,
                        alt,
                        state,
                        vgt,
                        gq,
                        bq,
                        read_depth,
                        ref_count,
                        alt_count,
                        vaf,
                    )
                ) 
            else:
                o.write(
                    "{}\t{}\t.\t{}\t{}\t.\t{}\t.\tGT:GQ:BQ:DP:AD:VAF\t{}:{}:{:0.1f}:{:0.0f}:{:0.0f},{:0.0f}:{:.2f}\n".format(
                        chrom,
                        pos,
                        ref,
                        alt,
                        state,
                        vgt,
                        gq,
                        bq,
                        read_depth,
                        ref_count,
                        alt_count,
                        vaf,
                    )
                )
    o.close()



def dump_snp_log(
    chrom_lst: List[str], 
    chrom2snp_log: Dict[str, List[int]], 
) -> None:

    row_names = [
        "num_snp",
        "num_low_gq_snp",
        "num_uncallable_snp",
        "num_ab_filtered_snp",
        "num_md_filtered_snp",
        "num_het_snp",
        "num_hetalt_snp",
        "num_homalt_snp",
    ]
    ncol = len(chrom_lst)
    nrow = len(row_names)
    dt = np.zeros((nrow, ncol))
    for i, chrom in enumerate(chrom_lst):
        for j, count in enumerate(chrom2snp_log[chrom]): 
            dt[j][i] = count
    
    o = open("himut.log", "w")
    col_lst = chrom_lst + ["total"] 
    o.write("{:30}{}\n".format("", "\t".join(col_lst)))
    for k in range(nrow):
        rsum =  str(int(np.sum(dt[k])))
        rlst = [str(int(r)) for r in dt[k].tolist()] + [rsum]
        o.write("{:30}{}\n".format(row_names[k], "\t".join(rlst)))
    o.close()