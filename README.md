## User Guide

ccs2snp is a basic germline mutation detection algorithm. ccs2snp leverages the base base accuracy and read length of Pacific Biosciences (PacBio) CCS reads to call germline SNPs.

## Usage

```sh
## use ccs2snp to call single nucleotide polymorphisms

ccs2snp call -i sample.primary_alignments.sorted.bam --ref ref.fasta -o sample.vcf  --region_list ref.target 
```

