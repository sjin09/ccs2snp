import time
import pysam
import pyfastx
import himut.util
import himut.bamlib
import himut.caller
import ccs2snp.gtlib
import ccs2snp.vcflib
import numpy as np
import multiprocessing as mp
from collections import defaultdict
from typing import Dict, List, Tuple



def init_allelecounts():
    rpos2allelecounts = defaultdict(lambda: np.zeros(6))
    rpos2allele2bq_lst = defaultdict(lambda: {0: [], 1: [], 2: [], 3: [], 4: [], 5: []})
    return rpos2allelecounts, rpos2allele2bq_lst


def update_allelecounts(
    ccs,
    rpos2allelecounts: Dict[int, np.ndarray],
    rpos2allele2bq_lst: Dict[int, Dict[int, List[int]]],
):

    tpos = ccs.tstart
    qpos = ccs.qstart
    for (state, ref, alt, ref_len, alt_len) in ccs.cstuple_lst:
        if state == 1:  # match
            for i, alt_base in enumerate(alt):
                epos = tpos + i
                bidx = himut.util.base2idx[alt_base]
                rpos2allelecounts[epos][bidx] += 1
                rpos2allele2bq_lst[epos][bidx].append(ccs.bq_int_lst[qpos + i])
        elif state == 2:  # sub
            bidx = himut.util.base2idx[alt]
            rpos2allelecounts[tpos][bidx] += 1
            rpos2allele2bq_lst[tpos][bidx].append(ccs.bq_int_lst[qpos])
        elif state == 3:  # insertion
            rpos2allelecounts[tpos][4] += 1
        elif state == 4:  # deletion
            for j in range(len(ref[1:])):
                rpos2allelecounts[tpos + j][5] += 1
        tpos += ref_len
        qpos += alt_len

    
def get_alt_counts(
    gt: str, 
    gt_state: str,
    read_depth: int,
    allelecounts: Dict[int, Dict[int, int]],
    allele2bq_lst: Dict[int, Dict[str, List[int]]],
):

    if gt_state == "hetalt":
        a1, a2 = list(gt)
        alt = ",".join([a1, a2])
        a1_idx = himut.util.base2idx[a1]
        a2_idx = himut.util.base2idx[a2]
        a1_count = allelecounts[a1_idx]
        a2_count = allelecounts[a2_idx]
        a1_bq = sum(allele2bq_lst[a1_idx])/float(a1_count)
        a2_bq = sum(allele2bq_lst[a2_idx])/float(a2_count)
        alt_bq = "{:0.1f},{:0.1f}".format(a1_bq, a2_bq)
        alt_count = "{:0.0f},{:0.0f}".format(a1_count, a2_count)
        alt_vaf = "{:.2f},{:.2f}".format(
            a1_count/float(read_depth), a2_count/float(read_depth)
        )
    elif gt_state == "het" or gt_state == "homalt":
        alt = gt[1]
        alt_count = allelecounts[himut.util.base2idx[alt]]
        alt_bq = sum(allele2bq_lst[himut.util.base2idx[alt]])/float(alt_count)
        alt_vaf = alt_count/float(read_depth)
    return alt, alt_bq, alt_vaf, alt_count



def get_germline_substitutions(
    chrom: str,
    bam_file: str,
    seq: str,
    chunkloci_lst: List[Tuple[str, int, int]],
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
    chrom2snp_lst: Dict[
        str, List[Tuple[str, int, str, str, int, int, int, int, str]]
    ],
):


    snp_lst = []
    himut.gtlib.init(germline_snv_prior)
    alignments = pysam.AlignmentFile(bam_file, "rb")
    for (chrom, chunk_start, chunk_end) in chunkloci_lst: 
        rpos2allelecounts, rpos2allele2bq_lst, = init_allelecounts()
        for i in alignments.fetch(chrom, chunk_start, chunk_end): # traverse the chromosome in chunks
            ccs = himut.bamlib.BAM(i)
            if not ccs.is_primary:
                continue
            if himut.caller.is_low_mapq(ccs.mapq, min_mapq):
                continue
            update_allelecounts(
                ccs, rpos2allelecounts, rpos2allele2bq_lst, 
            )

        for rpos in range(chunk_start, chunk_end):
            ref = seq[rpos]
            if ref not in himut.util.base_set:
                continue
            allelecounts = rpos2allelecounts[rpos]
            allele2bq_lst = rpos2allele2bq_lst[rpos]
            gt, gq, vgt, gt_state = ccs2snp.gtlib.get_germ_gt(ref, allele2bq_lst)
            if gt_state == "homref":
                continue

            tpos = rpos + 1
            _del_count, _ins_count, read_depth = himut.bamlib.get_read_depth(allelecounts)
            _, _, ref_count = himut.bamlib.get_ref_counts(ref, read_depth, allelecounts, allele2bq_lst)
            alt, alt_bq, alt_vaf, alt_count = get_alt_counts(gt, gt_state, read_depth, allelecounts, allele2bq_lst)
            if himut.caller.is_low_gq(gq, min_gq):
                snp_lst.append(
                    (
                        chrom,
                        tpos,
                        ref,
                        alt,
                        "LowGQ",
                        vgt,
                        gq,
                        alt_bq,
                        read_depth,
                        ref_count,
                        alt_count,
                        alt_vaf,
                    )
                )
                continue 

            if read_depth > md_threshold:
                snp_lst.append(
                    (
                        chrom,
                        tpos,
                        ref,
                        alt,
                        "HighDepth",
                        vgt,
                        gq,
                        alt_bq,
                        read_depth,
                        ref_count,
                        alt_count,
                        alt_vaf,
                    )
                )
                continue
            if gt_state == "het" and not (ref_count >= min_ref_count and alt_count >= min_alt_count):
                snp_lst.append(
                    (
                        chrom,
                        tpos,
                        ref,
                        alt,
                        "LowDepth",
                        vgt,
                        gq,
                        alt_bq,
                        read_depth,
                        ref_count,
                        alt_count,
                        alt_vaf,
                    )
                )
                continue

            counter = 0
            for j in alignments.fetch(chrom, tpos, tpos + 1):
                ccs = himut.bamlib.BAM(j)
                if not ccs.is_primary:
                    continue
                ccs.cs2subindel()
                tpos_lst = [tsbs[0] for tsbs in ccs.tsbs_lst]
                tpos_set = set(tpos_lst)
                if tpos not in tpos_set:
                    continue
                qpos = [qsbs[0] for qsbs in ccs.qsbs_lst][tpos_lst.index(tpos)]
                if himut.bamlib.is_mismatch_conflict(ccs, tpos, qpos, mismatch_window, max_mismatch_count):
                    counter += 1

            if alt_count == counter:
                snp_lst.append(
                    (
                        chrom,
                        tpos,
                        ref,
                        alt,
                        "MismatchConflict",
                        vgt,
                        gq,
                        alt_bq,
                        read_depth,
                        ref_count,
                        alt_count,
                        alt_vaf,
                    )
                )
                continue
            
            snp_lst.append(
                    (
                        chrom,
                        tpos,
                        ref,
                        alt,
                        "PASS",
                        vgt,
                        gq,
                        alt_bq,
                        read_depth,
                        ref_count,
                        alt_count,
                        alt_vaf,
                    )
                )
 
    chrom2snp_lst[chrom] = snp_lst
    alignments.close()


def call_germline_snps(
    bam_file: str,
    ref_file: str,
    region: str,
    region_list: str,
    min_bq: int,
    min_gq: int,
    min_mapq: int,
    min_ref_count: int,
    min_alt_count: int,
    mismatch_window: int,
    max_mismatch_count: int,
    germline_snv_prior: float,
    germline_indel_prior: float,
    threads: int,
    version: str,
    out_file: str,
) -> None:

    cpu_start = time.time() / 60
    ccs2snp.util.is_input_exist(
        bam_file,
        ref_file,
        out_file,
    )
    _, tname2tsize = himut.bamlib.get_tname2tsize(bam_file)
    chrom_lst, chrom2chunkloci_lst = himut.util.load_loci(region, region_list, tname2tsize)
    _, _, md_threshold = himut.bamlib.get_thresholds(bam_file, chrom_lst, tname2tsize)
    ccs2snp.util.is_input_corrupt(
        bam_file,
        ref_file,
        chrom_lst,
    )

    print("ccs2snp is calling substitutions with {} threads".format(threads))
    p = mp.Pool(threads)
    manager = mp.Manager()
    chrom2snp_lst = manager.dict()
    refseq = pyfastx.Fasta(ref_file)
    get_germline_substitutions_arg_lst = [
        (
            chrom,
            bam_file,
            str(refseq[chrom]),
            chrom2chunkloci_lst[chrom],
            min_bq,
            min_gq,
            min_mapq,
            md_threshold,
            min_ref_count,
            min_alt_count,
            mismatch_window,
            max_mismatch_count,
            germline_snv_prior,
            germline_indel_prior,
            chrom2snp_lst,
        )
        for chrom in chrom_lst
    ]
    p.starmap(
        get_germline_substitutions,
        get_germline_substitutions_arg_lst,
    )
    p.close()
    p.join()

    vcf_header = ccs2snp.vcflib.get_vcf_header(
        bam_file,
        ref_file,
        region,
        region_list,
        tname2tsize,
        min_bq, 
        min_gq,
        min_mapq,
        md_threshold,
        min_ref_count,
        min_alt_count,
        mismatch_window,
        max_mismatch_count,
        germline_snv_prior,
        germline_indel_prior,
        threads,
        version,
        out_file,
    )
    ccs2snp.vcflib.dump_snp(out_file, vcf_header, chrom_lst, chrom2snp_lst)

    print(
        "ccs2snp finished calling and returning germline SNPs with {} threads".format(
            threads
        )
    )
    cpu_end = time.time() / 60
    duration = cpu_end - cpu_start
    print(
        "ccs2snp germline SNP detection took {} minutes".format(
            duration
        )
    )
    himut.util.exit()
