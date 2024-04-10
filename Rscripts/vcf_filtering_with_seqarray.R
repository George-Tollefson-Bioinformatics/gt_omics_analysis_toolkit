library(gdsfmt)
library(SeqArray)

################################################################################
# Development Notes:
# this will be executable and also plugged into the shiny app (maybe loaded into shiny app via library)
#
# Decide whether to calculate missingness across sites after applying sample missingness filter (could start by getting mean missingness across unfiltered site and samples, and filter by whichever has higher missingness first to retain as much data as possible. 
# This may be question specific as well if we have gvcf with key mutations, we want to apply site missingness before sample missingness since we'd lose all samples even though they have our key mutations)
#
################################################################################

# Specify the path to your VCF file
vcf_file_path <- "/Users/george/scratch/prx07PRX07_biallelicsnps.maxsnpmissing_0.5.q30.sample_0.9.recode.vcf.gz"

# Specify the output path for the GDS file
gds_file_path <- "/Users/george/scratch/prx07_test_filtering_script.gds"

# Convert VCF to GDS format
seqVCF2GDS(vcf_file_path, gds_file_path)

# Open the newly created GDS file
gds <- seqOpen(gds_file_path) # this is a genotype file

# For efficiency, delete any fields which w know we will not sure - can change depending on analysis settings or tool used.
# delete "HM2", "HM3", "AA", "OR" and "DP"
# seqDelete(gds, info.var=c("HM2", "HM3", "AA", "OR"), fmt.var="DP")

### If user supplies a site missingness threshold
# filter open gds object in place. If changing filter, would need to close the gds object and re-open the gds object and filter, or store/write intermediate filtered gds objects
# Calculate the sample missingness for each variant
site_missingness <- seqMissing(gds, per.variant = TRUE) # FALSE runs missingness per sample - variants and samples is NA

# Identify variants with more than 90% missingness
variants_to_include <- which(site_missingness <= 0.9)
variant_ids_to_include <- seqGetData(gds, "variant.id")[variants_to_include]
seqSetFilter(gds, variant.id=variant_ids_to_include)

################################################################################
#### If user supplies a sample missingness threshold 
################################################################################

# filter open gds object in place. If changing filter, would need to close the gds object and re-open the gds object and filter, or store/write intermediate filtered gds objects
sample_missingness <- seqMissing(gds, per.variant = FALSE) # FALSE runs missingness per sample - variants and samples is NA

# Identify samples with more than 50% missingness
samples_to_include <- which(sample_missingness <= 0.5)

sample_ids_to_include <- seqGetData(gds, "sample.id")[samples_to_include]

################################################################################
# If user supplies a sample exclusion list:
################################################################################
# Filter out samples provided in sample_list.txt (newline separated sample ids)
sample_ids_to_include <- seqGetData(gds, "sample.id")[samples_to_include]
sample_ids_to_include <- sample_ids_to_include[-sample_exclusion_list]

seqSetFilter(gds, sample.id=sample_ids_to_include)

# since seqGetData dioesn't allow exclusion filtering, need to remove user defined sample list from our samples to include object
# If user does not specify a sample missingness filter but does supply sample exclusion list

################################################################################
# If user supplies a sample exclusion list:
################################################################################
# Filter out samples provided in sample_list.txt (newline separated sample ids)
sample_ids_to_include_no_filter <- seqGetData(gds, "sample.id")
sample_ids_to_include_no_filter <- sample_ids_to_include_no_filter[-sample_exclusion_list]

# Plot missingness histograms for sample and site missingness

# Calculate distribution of samples across sample missingess bins
# Create a data frame for plotting
missingness_df <- data.frame(Missingness = sample_missingness)

# Define bins for the histogram
bins <- seq(0, 1, by=0.1)

# Adjusting labels to correctly match the bins
# The last bin will be "90%-100%" hence no need to go beyond 90 in the sequence for labels
labels <- paste0(seq(0, 90, by=10), "%-", seq(10, 100, by=10), "%")

# Categorize missingness into bins
missingness_df$Bin <- cut(missingness_df$Missingness, bins, labels = labels, include.lowest = TRUE)

# Plotting the histogram with ggplot2
ggplot(missingness_df, aes(x = Bin)) +
  geom_bar(stat = "count", fill = "steelblue") +
  labs(title = "Distribution of Sample Missingness", x = "Missingness Range", y = "Number of Samples") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Rotate x-axis labels for better readability

# Save the plot
# User defined output directory
output_dir <- "/Users/george/scratch"

# Define plot title dynamically using user defined parameters
ggsave(paste0(output_dir,"\",sample_missingness_distribution.png"), width = 10, height = 6)

# Close the GDS file
seqClose(gds)

# Can write filtered table: 
# Inform the user
# cat("Filtered GDS file created:", filtered_gds_file_path, "\n")