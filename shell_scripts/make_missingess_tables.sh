#!/bin/bash
# Script to process VCF for missingness and prepare data for plotting
# Usage: script.sh <vcf_path> <output_path>
# vcf_path is full vcf path
# output path is the path where the output files will be saved

# Load necessary modules
module load vcftools

# Check if two arguments are provided
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <VCF_PATH> <OUTPUT_PATH> <PREFIX_TO_USE>"
    exit 1
fi

# Assign arguments to variables
vcf_path=$1
output_path=$2
prefix_to_use=$3

# Ensure output path ends with a slash
output_path=$(echo "$output_path" | sed 's:/*$::')/

# Assign variables for ease of use
vcf_to_use="$vcf_path"
output_with_prefix="${output_path}${prefix_to_use}"

# Calculate individual missingness
vcftools --gzvcf "$vcf_to_use" --missing-indv --out "$output_with_prefix"

# Prepare sample missingness table for R plotting
awk 'BEGIN { OFS="\t"; print "MissingnessThreshold", "SampleCount" }
     NR > 1 {
         miss = $5 * 100; # Convert fraction to percentage
         bin = int(miss / 10) * 10; # Find the bin (floor to nearest 10%)
         count[bin]++; # Increment the count for this bin
     }
     END {
         for (i = 0; i <= 90; i += 10) {
             print i "%-" i+10 "%", count[i]+0; # print count for each bin, adding 0 handles undefined bins
         }
         print "90%-100%", count[100]+0; # explicitly handle the 90-100% bin
     }' "${output_with_prefix}.imiss" > "${output_with_prefix}_sample_missingness_histogram_data.tsv"

# Calculate site missingness
vcftools --gzvcf "$vcf_to_use" --missing-site --out "$output_with_prefix"

awk 'BEGIN { OFS="\t"; print "MissingnessThreshold", "SiteCount" }
     NR > 1 {
         miss = $6 * 100; # Convert fraction to percentage
         bin = int(miss / 10) * 10; # find the bin (floor to nearest 10%)
         if (miss == 100) { # handle 100% missingness explicitly
             bin = 90;
         }
         count[bin]++; # increment the count for this bin
     }
     END {
         for (i = 0; i < 90; i += 10) {
             printf "%d%%-%d%%\t%d\n", i, i+10, count[i]+0;
         }
         printf "90%%-100%%\t%d\n", count[90]+0; # explicitly handle the 90-100% bin
     }' "${output_with_prefix}.lmiss" > "${output_with_prefix}_site_missingness_histogram_data.tsv"