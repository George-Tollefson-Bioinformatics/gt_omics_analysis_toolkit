# gt_omics_analysis_toolkit
A growing collection of custom executables for bioinformatics analysis

### Filter VCF:



Filtering a VCF file by include biallelic SNPs only, with minimum genotype quality, sample and snp missingness with optional sample_list filter specified by user.

```bash
Usage: bash shell_scripts/filter_vcf.sh <filename.vcf.gz> <max SNP missingness> <max sample missingness> <genotype quality> <output directory> <sample prefix> [sample list.txt]"
```

### Plot SNP and Sample Missingness:

First, make missingness tables. Run:
```bash
bash filter_vcf.sh <filename.vcf.gz> <max SNP missingness> <max sample missingness> <genotype quality> <output directory> <sample prefix> [sample list.txt]"
```

Then plot the sample and snp missingness histograms. Run:
```
module load r

Rscript Rscripts/plot_site_and_sample_missingness.R prefix <sample_prefix> <output_path> <sample_missingness_histogram_data.tsv> <snp_missingness_histogram_data.tsv>

```
