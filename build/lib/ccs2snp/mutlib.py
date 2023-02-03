import cyvcf2
import pyfastx
import himut.util
import himut.vcflib
import pandas as pd
# from plotnine import *
from collections import defaultdict
from typing import Dict, List, Tuple

purine = set(["A", "G"])
pyrimidine = set(["T", "C"])
sub_lst = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]
purine2pyrimidine = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}
sbs_lst = [
    "A[C>A]A",
    "A[C>A]C",
    "A[C>A]G",
    "A[C>A]T",
    "A[C>G]A",
    "A[C>G]C",
    "A[C>G]G",
    "A[C>G]T",
    "A[C>T]A",
    "A[C>T]C",
    "A[C>T]G",
    "A[C>T]T",
    "A[T>A]A",
    "A[T>A]C",
    "A[T>A]G",
    "A[T>A]T",
    "A[T>C]A",
    "A[T>C]C",
    "A[T>C]G",
    "A[T>C]T",
    "A[T>G]A",
    "A[T>G]C",
    "A[T>G]G",
    "A[T>G]T",
    "C[C>A]A",
    "C[C>A]C",
    "C[C>A]G",
    "C[C>A]T",
    "C[C>G]A",
    "C[C>G]C",
    "C[C>G]G",
    "C[C>G]T",
    "C[C>T]A",
    "C[C>T]C",
    "C[C>T]G",
    "C[C>T]T",
    "C[T>A]A",
    "C[T>A]C",
    "C[T>A]G",
    "C[T>A]T",
    "C[T>C]A",
    "C[T>C]C",
    "C[T>C]G",
    "C[T>C]T",
    "C[T>G]A",
    "C[T>G]C",
    "C[T>G]G",
    "C[T>G]T",
    "G[C>A]A",
    "G[C>A]C",
    "G[C>A]G",
    "G[C>A]T",
    "G[C>G]A",
    "G[C>G]C",
    "G[C>G]G",
    "G[C>G]T",
    "G[C>T]A",
    "G[C>T]C",
    "G[C>T]G",
    "G[C>T]T",
    "G[T>A]A",
    "G[T>A]C",
    "G[T>A]G",
    "G[T>A]T",
    "G[T>C]A",
    "G[T>C]C",
    "G[T>C]G",
    "G[T>C]T",
    "G[T>G]A",
    "G[T>G]C",
    "G[T>G]G",
    "G[T>G]T",
    "T[C>A]A",
    "T[C>A]C",
    "T[C>A]G",
    "T[C>A]T",
    "T[C>G]A",
    "T[C>G]C",
    "T[C>G]G",
    "T[C>G]T",
    "T[C>T]A",
    "T[C>T]C",
    "T[C>T]G",
    "T[C>T]T",
    "T[T>A]A",
    "T[T>A]C",
    "T[T>A]G",
    "T[T>A]T",
    "T[T>C]A",
    "T[T>C]C",
    "T[T>C]G",
    "T[T>C]T",
    "T[T>G]A",
    "T[T>G]C",
    "T[T>G]G",
    "T[T>G]T",
]
tri_lst = [
    "ACA",
    "ACC",
    "ACG",
    "ACT",
    "ATA",
    "ATC",
    "ATG",
    "ATT",
    "CCA",
    "CCC",
    "CCG",
    "CCT",
    "CTA",
    "CTC",
    "CTG",
    "CTT",
    "GCA",
    "GCC",
    "GCG",
    "GCT",
    "GTA",
    "GTC",
    "GTG",
    "GTT",
    "TCA",
    "TCC",
    "TCG",
    "TCT",
    "TTA",
    "TTC",
    "TTG",
    "TTT",
]
sbs2sub = {
    "A[C>A]A": "C>A",
    "A[C>A]T": "C>A",
    "A[C>A]G": "C>A",
    "A[C>A]C": "C>A",
    "T[C>A]A": "C>A",
    "T[C>A]T": "C>A",
    "T[C>A]G": "C>A",
    "T[C>A]C": "C>A",
    "G[C>A]A": "C>A",
    "G[C>A]T": "C>A",
    "G[C>A]G": "C>A",
    "G[C>A]C": "C>A",
    "C[C>A]A": "C>A",
    "C[C>A]T": "C>A",
    "C[C>A]G": "C>A",
    "C[C>A]C": "C>A",
    "A[C>G]A": "C>G",
    "A[C>G]T": "C>G",
    "A[C>G]G": "C>G",
    "A[C>G]C": "C>G",
    "T[C>G]A": "C>G",
    "T[C>G]T": "C>G",
    "T[C>G]G": "C>G",
    "T[C>G]C": "C>G",
    "G[C>G]A": "C>G",
    "G[C>G]T": "C>G",
    "G[C>G]G": "C>G",
    "G[C>G]C": "C>G",
    "C[C>G]A": "C>G",
    "C[C>G]T": "C>G",
    "C[C>G]G": "C>G",
    "C[C>G]C": "C>G",
    "A[C>T]A": "C>T",
    "A[C>T]T": "C>T",
    "A[C>T]G": "C>T",
    "A[C>T]C": "C>T",
    "T[C>T]A": "C>T",
    "T[C>T]T": "C>T",
    "T[C>T]G": "C>T",
    "T[C>T]C": "C>T",
    "G[C>T]A": "C>T",
    "G[C>T]T": "C>T",
    "G[C>T]G": "C>T",
    "G[C>T]C": "C>T",
    "C[C>T]A": "C>T",
    "C[C>T]T": "C>T",
    "C[C>T]G": "C>T",
    "C[C>T]C": "C>T",
    "A[T>A]A": "T>A",
    "A[T>A]T": "T>A",
    "A[T>A]G": "T>A",
    "A[T>A]C": "T>A",
    "T[T>A]A": "T>A",
    "T[T>A]T": "T>A",
    "T[T>A]G": "T>A",
    "T[T>A]C": "T>A",
    "G[T>A]A": "T>A",
    "G[T>A]T": "T>A",
    "G[T>A]G": "T>A",
    "G[T>A]C": "T>A",
    "C[T>A]A": "T>A",
    "C[T>A]T": "T>A",
    "C[T>A]G": "T>A",
    "C[T>A]C": "T>A",
    "A[T>C]A": "T>C",
    "A[T>C]T": "T>C",
    "A[T>C]G": "T>C",
    "A[T>C]C": "T>C",
    "T[T>C]A": "T>C",
    "T[T>C]T": "T>C",
    "T[T>C]G": "T>C",
    "T[T>C]C": "T>C",
    "G[T>C]A": "T>C",
    "G[T>C]T": "T>C",
    "G[T>C]G": "T>C",
    "G[T>C]C": "T>C",
    "C[T>C]A": "T>C",
    "C[T>C]T": "T>C",
    "C[T>C]G": "T>C",
    "C[T>C]C": "T>C",
    "A[T>G]A": "T>G",
    "A[T>G]T": "T>G",
    "A[T>G]G": "T>G",
    "A[T>G]C": "T>G",
    "T[T>G]A": "T>G",
    "T[T>G]T": "T>G",
    "T[T>G]G": "T>G",
    "T[T>G]C": "T>G",
    "G[T>G]A": "T>G",
    "G[T>G]T": "T>G",
    "G[T>G]G": "T>G",
    "G[T>G]C": "T>G",
    "C[T>G]A": "T>G",
    "C[T>G]T": "T>G",
    "C[T>G]G": "T>G",
    "C[T>G]C": "T>G",
}
sbs2tri = {
    "A[C>A]A": "ACA",
    "A[C>A]T": "ACT",
    "A[C>A]G": "ACG",
    "A[C>A]C": "ACC",
    "T[C>A]A": "TCA",
    "T[C>A]T": "TCT",
    "T[C>A]G": "TCG",
    "T[C>A]C": "TCC",
    "G[C>A]A": "GCA",
    "G[C>A]T": "GCT",
    "G[C>A]G": "GCG",
    "G[C>A]C": "GCC",
    "C[C>A]A": "CCA",
    "C[C>A]T": "CCT",
    "C[C>A]G": "CCG",
    "C[C>A]C": "CCC",
    "A[C>G]A": "ACA",
    "A[C>G]T": "ACT",
    "A[C>G]G": "ACG",
    "A[C>G]C": "ACC",
    "T[C>G]A": "TCA",
    "T[C>G]T": "TCT",
    "T[C>G]G": "TCG",
    "T[C>G]C": "TCC",
    "G[C>G]A": "GCA",
    "G[C>G]T": "GCT",
    "G[C>G]G": "GCG",
    "G[C>G]C": "GCC",
    "C[C>G]A": "CCA",
    "C[C>G]T": "CCT",
    "C[C>G]G": "CCG",
    "C[C>G]C": "CCC",
    "A[C>T]A": "ACA",
    "A[C>T]T": "ACT",
    "A[C>T]G": "ACG",
    "A[C>T]C": "ACC",
    "T[C>T]A": "TCA",
    "T[C>T]T": "TCT",
    "T[C>T]G": "TCG",
    "T[C>T]C": "TCC",
    "G[C>T]A": "GCA",
    "G[C>T]T": "GCT",
    "G[C>T]G": "GCG",
    "G[C>T]C": "GCC",
    "C[C>T]A": "CCA",
    "C[C>T]T": "CCT",
    "C[C>T]G": "CCG",
    "C[C>T]C": "CCC",
    "A[T>A]A": "ATA",
    "A[T>A]T": "ATT",
    "A[T>A]G": "ATG",
    "A[T>A]C": "ATC",
    "T[T>A]A": "TTA",
    "T[T>A]T": "TTT",
    "T[T>A]G": "TTG",
    "T[T>A]C": "TTC",
    "G[T>A]A": "GTA",
    "G[T>A]T": "GTT",
    "G[T>A]G": "GTG",
    "G[T>A]C": "GTC",
    "C[T>A]A": "CTA",
    "C[T>A]T": "CTT",
    "C[T>A]G": "CTG",
    "C[T>A]C": "CTC",
    "A[T>C]A": "ATA",
    "A[T>C]T": "ATT",
    "A[T>C]G": "ATG",
    "A[T>C]C": "ATC",
    "T[T>C]A": "TTA",
    "T[T>C]T": "TTT",
    "T[T>C]G": "TTG",
    "T[T>C]C": "TTC",
    "G[T>C]A": "GTA",
    "G[T>C]T": "GTT",
    "G[T>C]G": "GTG",
    "G[T>C]C": "GTC",
    "C[T>C]A": "CTA",
    "C[T>C]T": "CTT",
    "C[T>C]G": "CTG",
    "C[T>C]C": "CTC",
    "A[T>G]A": "ATA",
    "A[T>G]T": "ATT",
    "A[T>G]G": "ATG",
    "A[T>G]C": "ATC",
    "T[T>G]A": "TTA",
    "T[T>G]T": "TTT",
    "T[T>G]G": "TTG",
    "T[T>G]C": "TTC",
    "G[T>G]A": "GTA",
    "G[T>G]T": "GTT",
    "G[T>G]G": "GTG",
    "G[T>G]C": "GTC",
    "C[T>G]A": "CTA",
    "C[T>G]T": "CTT",
    "C[T>G]G": "CTG",
    "C[T>G]C": "CTC",
}


def get_sbs48(
    chrom: str,
    pos: str,
    ref: str,
    alt: str,
    refseq,
) -> str:
    if ref in purine:
        upstream = purine2pyrimidine.get(refseq[chrom][pos + 1], "N")
        downstream = purine2pyrimidine.get(refseq[chrom][pos - 1], "N")
        sbs96 = "{}[{}>{}]{}".format(
            upstream,
            purine2pyrimidine.get(ref, "N"),
            purine2pyrimidine.get(alt, "N"),
            downstream,
        )
    else:
        upstream = refseq[chrom][pos - 1]
        downstream = refseq[chrom][pos + 1]
        sbs96 = "{}[{}>{}]{}".format(upstream, ref, alt, downstream)
    return sbs96


def load_sbs48_counts(
    vcf_file: str, 
    ref_file: str, 
    chrom_lst: List[str]
) -> Dict[str, int]:

    tname_lst = []
    tname2sbs2counts = defaultdict()
    refseq = pyfastx.Fasta(ref_file)
    if vcf_file.endswith(".vcf"):
        for line in open(vcf_file).readlines():
            if line.startswith("##"):
                if line.startswith("##contig"):
                    tname = line.strip().replace("##contig=<ID=", "").split(",")[0]
                    tname_lst.append(tname)
                continue
            elif line.startswith("#CHROM"):
                for tname in tname_lst:
                    tname2sbs2counts[tname] = defaultdict(lambda: 0)
                continue
            v = himut.vcflib.VCF(line)
            if v.is_snp and v.is_pass:
                sbs96 = get_sbs96(v.chrom, int(v.pos) - 1, v.ref, v.alt, refseq)
                tname2sbs2counts[v.chrom][sbs96] += 1
    elif vcf_file.endswith(".vcf.bgz"):
        for line in cyvcf2.VCF(vcf_file).raw_header.split("\n"):
            if line.startswith("##"):
                if line.startswith("##contig"):
                    tname = line.replace("##contig=<ID=", "").split(",")[0]
                    tname_lst.append(tname)
                continue
            elif line.startswith("#CHROM"):
                for tname in tname_lst:
                    tname2sbs2counts[tname] = defaultdict(lambda: 0)
        for i in cyvcf2.VCF(vcf_file):
            v = himut.vcflib.VCF(str(i))
            if v.is_snp and v.is_pass:
                sbs96 = get_sbs96(v.chrom, int(v.pos) - 1, v.ref, v.alt, refseq)
                tname2sbs2counts[v.chrom][sbs96] += 1

    sbs2count = defaultdict(lambda: 0)
    for chrom in chrom_lst:
        for sbs, count in tname2sbs2counts[chrom].items():
            sbs2count[sbs] += count
    return sbs2count


def dump_sbs48_counts(
    vcf_file: str, ref_file: str, chrom: str, chrom_fofn: str, out_file: str
) -> None:

    chrom2sbs2counts = load_sbs96_counts(vcf_file, ref_file)
    chrom_lst = himut.util.load_chrom(chrom, chrom_fofn)
    sbs2counts = defaultdict(lambda: 0)
    if len(chrom_lst) == 0:
        for chrom in chrom2sbs2counts:
            for sbs, counts in chrom2sbs2counts[chrom].items():
                sbs2counts[sbs] += counts
    else:
        for chrom in chrom_lst:
            for sbs, counts in chrom2sbs2counts[chrom].items():
                sbs2counts[sbs] += counts

    o = open(out_file, "w")
    o.write("{}\t{}\t{}\t{}\n".format("sub", "tri", "sbs96", "counts"))
    for sbs in sbs_lst:
        o.write(
            "{}\t{}\t{}\t{}\n".format(sbs2sub[sbs], sbs2tri[sbs], sbs, sbs2counts[sbs])
        )
    o.close()


def dump_sbs48_plt(infile: str, sample: str, outfile: str) -> None:

    if sample is None:
        print(
            "Sample cannot be None\nBAM or VCF file might be missing sample information"
        )
        himut.util.exit()

    df = pd.read_csv(infile, sep="\t")
    plot = (
        ggplot(df, aes(x="tri", y="counts", fill="sub"))
        + geom_bar(stat="identity")
        + theme_bw()
        + facet_grid(". ~ sub", scales="free")
        + scale_fill_manual(
            values=("#98D7EC", "#212121", "#FF003A", "#A6A6A6", "#83A603", "#F5ABCC")
        )
        + labs(x="\nTrinucleotide Context\n", y="\nCounts\n")
        + ggtitle(sample)
        + theme(
            legend_title=element_blank(),
            axis_text_x=element_text(
                family="monospace", size=10, angle=90, ha="center"
            ),
            text=element_text(size=12),
        )
    )
    plot.save(outfile, width=22, height=10)

