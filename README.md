# gt_omics_analysis_toolkit
A growing collection of custom executables for bioinformatics analysis

# Filter VCF:

### shell_scripts/filter_vcf.sh

Filtering a VCF file by include biallelic SNPs only, with minimum genotype quality, sample and snp missingness with optional sample_list filter specified by user.

```bash
Usage: bash filter_vcf.sh <filename.vcf.gz> <max SNP missingness> <max sample missingness> <genotype quality> <output directory> <sample prefix> [sample list.txt]"
```

# Plot SNP and Sample Missingness:

### shell_scripts/filter_vcf.sh and Rscripts/plot_site_and_sample_missingness.R

First run:
```bash
Usage: bash filter_vcf.sh <filename.vcf.gz> <max SNP missingness> <max sample missingness> <genotype quality> <output directory> <sample prefix> [sample list.txt]"
```

Then run:
```
module load r

Usage: Rscript Rscripts/plot_site_and_sample_missingness.R prefix <sample_prefix> <output_path> <sample_missingness_histogram_data.tsv> <snp_missingness_histogram_data.tsv>

```
