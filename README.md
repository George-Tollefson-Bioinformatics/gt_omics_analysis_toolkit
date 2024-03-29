# gt_omics_analysis_toolkit
A growing collection of custom executables for bioinformatics analysis

# Shell Scripts:

### filter_vcf.sh

Filtering a VCF file by include biallelic SNPs only, with minimum genotype quality, sample and snp missingness with optional sample_list filter specified by user.

```bash
Usage: bash filter_vcf.sh <filename.vcf.gz> <max SNP missingness> <max sample missingness> <genotype quality> <output directory> <sample prefix> [sample list.txt]"
```
