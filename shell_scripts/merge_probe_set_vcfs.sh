#!/bin/bash
#SBATCH --job-name=merge_IBC_DR2_vcfs
#SBATCH --time=24:00:00 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --array=0-4
#SBATCH --nodes=1
#SBATCH --mem=32G

#####################################################################################################################
#
# Cross Probe Set VCF Merger - Concise, fast, and clean cross probe set VCF merging run via job array.
# George Tollefson April 2024
#
# Utilizes multithreading for bcftools I/O compression and decompression wherever possible.
#
# This script merges two VCF files from two different collections (DR2 and IBC) for each year of collection.
# The script first removes samples which are not in the intersection of the two VCFs, 
# then sorts, indexes, and merges the VCF files.
# The merged VCF file is then normalized to collapse any overlapping regions.
#
# It was originally written to merge our probe set VCFs for the DR2 and IBC collections for the Ugandan k13 2016-2022 analysis.
#####################################################################################################################
module load bcftools
module load tabix

# Define VCF files to merge. Merging two probe set vcfs for each collection year.
PRX00_DR2_vcf="/DR2/PRX00_current/variants/variants.vcf.gz"
PRX00_IBC_vcf="/variant_calling/IBC/PRX00_smk/variants.vcf.gz"

PRX02_DR2_vcf="/DR2/PRX02_current/variants/variants.vcf.gz"
PRX02_IBC_vcf="/variant_calling/IBC/PRX02_smk/variants.vcf.gz"

PRX05_DR2_vcf="/DR2/PRX05_current/variants/variants.vcf.gz"
PRX05_IBC_vcf="/variant_calling/IBC/PRX05/variants.vcf.gz"

PRX06_DR2_vcf="/DR2/PRX06_current/variants/variants.vcf.gz"
PRX06_IBC_vcf="/variant_calling/IBC/PRX06_smk/variants.vcf.gz"

PRX07_DR2_vcf="/DR2/PRX07_current/variants/variants.vcf.gz"
PRX07_IBC_vcf="/variant_calling/IBC/PRX07/variants.vcf.gz"

# Define array lists
DR2_vcf_list=($PRX00_DR2_vcf $PRX02_DR2_vcf $PRX05_DR2_vcf $PRX06_DR2_vcf $PRX07_DR2_vcf)
IBC_vcf_list=($PRX00_IBC_vcf $PRX02_IBC_vcf $PRX05_IBC_vcf $PRX06_IBC_vcf $PRX07_IBC_vcf)
prefix_list=("PRX00" "PRX02" "PRX05" "PRX06" "PRX07")

# Processing logic
vcf1=${DR2_vcf_list[$SLURM_ARRAY_TASK_ID]}
vcf2=${IBC_vcf_list[$SLURM_ARRAY_TASK_ID]}
prefix=${prefix_list[$SLURM_ARRAY_TASK_ID]}

# Define base directories
output_directory="/variant_calling/merged_vcfs/DR2_IBC"
temp_directory="${output_directory}/${prefix}_temp" # careful to make unqiue temp dir for each array task as they will not complete at same time and temp dir removal upon task completion will disrupt other tasks

# Create output directory if it doesn't exist
if [ ! -d "$output_directory" ]; then
    echo "Output directory does not exist, creating it..."
    mkdir -p "$output_directory"
fi

# Create temporary directory if it doesn't exist
mkdir -p $temp_directory

# Define temp and final output files
output_vcf_list=("PRX00_DR2_IBC_merged.vcf.gz" "PRX02_DR2_IBC_merged.vcf.gz" "PRX05_DR2_IBC_merged.vcf.gz" "PRX06_DR2_IBC_merged.vcf.gz" "PRX07_DR2_IBC_merged.vcf.gz")
norm_output_vcf_list=("PRX00_DR2_IBC_merged_norm.vcf.gz" "PRX02_DR2_IBC_merged_norm.vcf.gz" "PRX05_DR2_IBC_merged_norm.vcf.gz" "PRX06_DR2_IBC_merged_norm.vcf.gz" "PRX07_DR2_IBC_merged_norm.vcf.gz")

sample_intersect_vcf1="$temp_directory/${prefix}_$(basename $vcf1 .vcf.gz)_DR2_sorted.vcf.gz"
sample_intersect_vcf2="$temp_directory/${prefix}_$(basename $vcf2 .vcf.gz)_IBC_sorted.vcf.gz"

sorted_vcf1="$temp_directory/${prefix}_$(basename $vcf1 .vcf.gz)_DR2_sorted.vcf.gz"
sorted_vcf2="$temp_directory/${prefix}_$(basename $vcf2 .vcf.gz)_IBC_sorted.vcf.gz"

merged_vcf="$temp_directory/${output_vcf_list[$SLURM_ARRAY_TASK_ID]}"
sorted_merged_vcf="$temp_directory/${prefix}_$(basename $vcf2 .vcf.gz)_merged_sorted.vcf.gz"
norm_merged_vcf="$output_directory/${norm_output_vcf_list[$SLURM_ARRAY_TASK_ID]}"

# Remove samples which are not in the intersection of the two VCFs
echo "Removing samples not in the intersection of the two VCFs: $vcf1 and $vcf2"

# Extract sample names from each file
bcftools query -l $vcf1 > $temp_directory/${prefix}_samples1.txt
bcftools query -l $vcf2 > $temp_directory/${prefix}_samples2.txt


# Find common samples
grep -Fx -f $temp_directory/${prefix}_samples1.txt $temp_directory/${prefix}_samples2.txt > $temp_directory/${prefix}_common_samples.txt

# Filter both VCF files to include only common samples
bcftools view $vcf1 -S $temp_directory/${prefix}_common_samples.txt -Oz -o $sample_intersect_vcf1 --threads ${SLURM_CPUS_PER_TASK}
bcftools view $vcf2 -S $temp_directory/${prefix}_common_samples.txt -Oz -o $sample_intersect_vcf2 --threads ${SLURM_CPUS_PER_TASK}

# Index the filtered files
bcftools index $sample_intersect_vcf1
bcftools index $sample_intersect_vcf2

# Sort, index, and merge VCF files
echo "Sorting VCF files: $vcf1 and $vcf2"
bcftools sort $sample_intersect_vcf1 -Oz -o $sorted_vcf1
bcftools sort $sample_intersect_vcf2 -Oz -o $sorted_vcf2

echo "Indexing sorted VCF files: $sorted_vcf1 and $sorted_vcf2"
bcftools index $sorted_vcf1
bcftools index $sorted_vcf2

# Use concat not merge and follow with normalization to merge any variants regions whihc overlap in both starting vcf files
echo "Concatenating sorted VCF files: $sorted_vcf1 and $sorted_vcf2"
bcftools concat $sorted_vcf1 $sorted_vcf2 -a -Oz -o $merged_vcf --threads ${SLURM_CPUS_PER_TASK}
echo "Concatenation completed: $merged_vcf"

echo "Indexing concatenated VCF"
bcftools index $merged_vcf

# Sort again in case of any weird overalpping regions being concat out of order
echo "Sorting merged VCF"
bcftools sort $merged_vcf -Oz -o $sorted_merged_vcf

# Normalize to collapse any overlapping regions
echo "Normalizing merged VCF"
bcftools norm -m-any $sorted_merged_vcf -Oz -o $norm_merged_vcf --threads ${SLURM_CPUS_PER_TASK}
bcftools index $norm_merged_vcf

# Cleanup and complete
echo "Cleaning temp files..."
rm -r $temp_directory
echo "Merge completed: ${norm_merged_vcf} is ready."
