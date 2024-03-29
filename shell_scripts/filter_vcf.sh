#!/bin/bash
# Script for filtering one VCF file at a time to include biallelic SNPs only, with minimum genotype quality, sample and snp missingness with optional sample_list filter specified by user.
# Written to be run on Oscar cluster at Brown where vcftools bcftools tabix are all installed
# example run command for testing: 
# bash /nfs/jbailey5/baileyweb/gtollefs/toolshed/gt_omics_analysis_toolkit/shell_scripts/filter_PRX_collections.sh /nfs/jbailey5/baileyweb/gtollefs/uganda_k13_selection/variant_calling/IBC/PRX07/variants.vcf.gz 0.5 0.9 30 /nfs/jbailey5/baileyweb/gtollefs/uganda_k13_selection/variant_calling/filtering/PRX07/ PRX07 /nfs/jbailey5/baileyweb/gtollefs/uganda_k13_selection/variant_calling/filtering/scripts/prx07_control_samples_to_remove.txt

module load vcftools bcftools tabix

# Check if the correct number of arguments were provided
if [ "$#" -lt 5 ] || [ "$#" -gt 7 ]; then
    echo "Usage: bash $0 <filename.vcf.gz> <max SNP missingness> <max sample missingness> <genotype quality> <output directory> <sample prefix> [sample list.txt]"
    echo "Example: $0 /your/path/filename.vcf.gz 0.5 0.5 30 /your/output/directory/ sample_prefix [sample_list.txt]"
    exit 1
fi

# Assign arguments to variables with default values for missingness and quality
filename=$1
max_snp_missingness=${2:-50}
max_sample_missingness=${3:-50}
genotype_quality=${4:-30}
output_dir=$5
sample_prefix=$6
sample_list=${7-}

# Ensure output directory ends with a slash
output_dir=$(echo "$output_dir" | sed 's:/*$::')/

temp_dir=$(mktemp -d)

echo "Using temporary directory $temp_dir for intermediate files"

# Echo inputs for verification
echo "Processing file: $filename"
echo "Max SNP Missingness: $max_snp_missingness"
echo "Max Sample Missingness: $max_sample_missingness"
echo "Genotype Quality: $genotype_quality"
echo "Output Directory: $output_dir"
echo "Sample Prefix: $sample_prefix"
if [[ -n "$sample_list" ]]; then
    echo "Sample List: $sample_list"
fi

# dev assignments:
# filename="/nfs/jbailey5/baileyweb/gtollefs/uganda_k13_selection/variant_calling/IBC/PRX07/variants.vcf.gz"
# max_snp_missingness="0.5"
# max_sample_missingness="0.5"
# output_dir="/nfs/jbailey5/baileyweb/gtollefs/uganda_k13_selection/variant_calling/filtering/PRX07/"
# sample_prefix="PRX07"

# filter out indels
echo "############################################"
echo "Filtering out indels."
echo "############################################"
vcftools --gzvcf $filename --remove-indels --recode --recode-INFO-all --out ${temp_dir}${sample_prefix}_snps_only
# VCFtools is deprecated (no longer supported) in favor of BCFTools but not sure what equivalent BCFTools filtering commands are

# filter for biallelic snps only
echo "############################################"
echo "Filtering to keep biallelic snps only."
echo "############################################"
vcftools --vcf ${temp_dir}${sample_prefix}_snps_only.recode.vcf --min-alleles 2 --max-alleles 2 --out ${temp_dir}${sample_prefix}_biallelicsnps --recode

# snp missingness and minQ (You can use different thresholds based on your data and research questions
echo "############################################"
echo "Filtering to keep snps with max snp missingness of $max_snp_missingness and minimum genotype quality of $genotype_quality."
echo "############################################"
# vcftools --vcf  ${temp_dir}${sample_prefix}_biallelicsnps.recode.vcf --max-missing $max_snp_missingness  --minQ 30 --recode --recode-INFO-all --out ${temp_dir}${sample_prefix}_biallelicsnps.maxsnpmissing_${max_snp_missingness}.q30
vcftools --vcf  ${temp_dir}${sample_prefix}_biallelicsnps.recode.vcf --max-missing $max_snp_missingness --minQ $genotype_quality --recode --recode-INFO-all --out ${temp_dir}${sample_prefix}_biallelicsnps.maxsnpmissing_${max_snp_missingness}.q${genotype_quality}

# filter for sample misingness
echo "############################################"
echo "Filtering to keep snps with max sample missingness of $max_sample_missingness."
echo "############################################"
# calculate Sample missingness and produce total_missing report
vcftools --vcf ${temp_dir}${sample_prefix}_biallelicsnps.maxsnpmissing_${max_snp_missingness}.q${genotype_quality}.recode.vcf --missing-indv
awk '!/IN/' out.imiss | cut -f5 > ${temp_dir}${sample_prefix}_biallelicsnps_recode.maxsnpmissing_${max_snp_missingness}.q${genotype_quality}_totalmissing.txt

# 50% missingness based on the data and check how many samples remaining 
# older version without string interpolation # awk '$5 > 0.5' out.imiss | cut -f1 > ${filename}_biallelicsnps_recode.maxsnpmissing_${max_snp_missingness}.q${genotype_quality}_low50DP.indv
awk -v max_miss="$max_sample_missingness" '$5 > max_miss' out.imiss | cut -f1 > ${temp_dir}${sample_prefix}_biallelicsnps_recode.maxsnpmissing_${max_snp_missingness}.q${genotype_quality}.${max_sample_missingness}_indv
vcftools --vcf ${temp_dir}${sample_prefix}_biallelicsnps.maxsnpmissing_${max_snp_missingness}.q${genotype_quality}.recode.vcf --remove ${temp_dir}${sample_prefix}_biallelicsnps_recode.maxsnpmissing_${max_snp_missingness}.q${genotype_quality}.${max_sample_missingness}_indv --recode --recode-INFO-all --out ${temp_dir}${sample_prefix}_biallelicsnps.maxsnpmissing_${max_snp_missingness}.q${genotype_quality}.sample_${max_sample_missingness}

# if sample list is provided, use it to filter samples in vcf
echo "############################################"
echo "Filtering to keep remove samples in ${sample_list} if provided."
echo "############################################"
if [[ -n "$sample_list" ]]; then
    # use newline-separated sample list to filter samples with bcftools
    bcftools view -S $sample_list --force-samples ${temp_dir}${sample_prefix}_biallelicsnps.maxsnpmissing_${max_snp_missingness}.q${genotype_quality}.sample_${max_sample_missingness}.recode.vcf > ${temp_dir}${sample_prefix}.biallelicsnps.maxsnpmissing_${max_snp_missingness}.q${genotype_quality}.sample_${max_sample_missingness}.recode.sample_list_filtered.vcf
    echo "Filtered VCF by provided sample list: $sample_list"

    # zip
    echo "############################################"
    echo "gzipping vcf with bgzip."
    echo "############################################"
    bgzip -c ${temp_dir}${sample_prefix}.biallelicsnps.maxsnpmissing_${max_snp_missingness}.q${genotype_quality}.sample_${max_sample_missingness}.recode.sample_list_filtered.vcf > ${output_dir}${sample_prefix}.biallelicsnps.maxsnpmissing_${max_snp_missingness}.q${genotype_quality}.sample_${max_sample_missingness}.recode.sample_list_filtered.vcf.gz
    echo "Output file: ${output_dir}${sample_prefix}.biallelicsnps.maxsnpmissing_${max_snp_missingness}.q${genotype_quality}.sample_${max_sample_missingness}.recode.sample_list_filtered.vcf.gz"

    # index
    echo "############################################"
    echo "indexing vcf.gz with bcftools and tabix to produce .csi and .tbi index files."
    echo "############################################"
    bcftools index -c ${output_dir}${sample_prefix}.biallelicsnps.maxsnpmissing_${max_snp_missingness}.q${genotype_quality}.sample_${max_sample_missingness}.recode.sample_list_filtered.vcf.gz
    tabix -p vcf ${output_dir}${sample_prefix}.biallelicsnps.maxsnpmissing_${max_snp_missingness}.q${genotype_quality}.sample_${max_sample_missingness}.recode.sample_list_filtered.vcf.gz

fi

# zip
echo "############################################"
echo "gzipping vcf with bgzip."
echo "############################################"
bgzip -c ${temp_dir}${sample_prefix}_biallelicsnps.maxsnpmissing_${max_snp_missingness}.q${genotype_quality}.sample_${max_sample_missingness}.recode.vcf > ${output_dir}${sample_prefix}_biallelicsnps.maxsnpmissing_${max_snp_missingness}.q${genotype_quality}.sample_${max_sample_missingness}.recode.vcf.gz

# index
echo "############################################"
echo "indexing vcf.gz with bcftools and tabix to produce .csi and .tbi index files."
echo "############################################"
bcftools index -c ${output_dir}${sample_prefix}_biallelicsnps.maxsnpmissing_${max_snp_missingness}.q${genotype_quality}.sample_${max_sample_missingness}.recode.vcf.gz
tabix -p vcf ${output_dir}${sample_prefix}_biallelicsnps.maxsnpmissing_${max_snp_missingness}.q${genotype_quality}.sample_${max_sample_missingness}.recode.vcf.gz

# clean temp directory
rm -r "$temp_dir"
echo "Temporary files cleaned up."

# move logs
echo "############################################"
echo "Moving all logs to one log directory within the specified output directory."
echo "############################################"
mkdir -p ${output_dir}logs
mv *.log ${output_dir}logs