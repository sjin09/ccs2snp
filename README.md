## User Guide

ccs2snp is a basic germline mutation detection algorithm that uses PacBio CCS reads to call germline SNPs.

## Usage

```sh
## use ccs2snp to call germline single nucleotide polymorphisms

ccs2snp -i sample.primary_alignments.sorted.bam --ref ref.fasta --region_list ref.target -o sample.vcf  
```

