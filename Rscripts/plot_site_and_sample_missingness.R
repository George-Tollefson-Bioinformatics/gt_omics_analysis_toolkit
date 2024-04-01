#!/usr/bin/env Rscript

# R script to report sample and site missingness

# Load necessary libraries
library(ggplot2)

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if the correct number of arguments is provided
if(length(args) != 4) {
  stop("Four arguments are required: sample_prefix, output_path, sample_missingness_file, and site_missingness_file", call. = FALSE)
}

# Assign arguments to variables
sample_prefix <- args[1]
output_path <- args[2]
sample_missingness_file <- args[3]
site_missingness_file <- args[4]

# Check if filenames contain the expected words
if(!grepl("sample", sample_missingness_file)) {
  stop("The sample missingness file name does not contain the word 'sample'", call. = FALSE)
}

if(!grepl("site", site_missingness_file)) {
  stop("The site missingness file name does not contain the word 'site'", call. = FALSE)
}

# Load data
sample_missingness_data <- read.table(sample_missingness_file, header = TRUE, sep="\t")
site_missingness_data <- read.table(site_missingness_file, header = TRUE, sep="\t")

# Plotting Sample Missingness
sample_plot <- ggplot(sample_missingness_data, aes(x=MissingnessThreshold, y=SampleCount)) +
  geom_bar(stat="identity", fill="skyblue") +
  geom_text(aes(label=SampleCount), vjust=-0.5, color="black") +
  theme_minimal() +
  labs(title="Histogram of Sample Missingness",
       x="Missingness Threshold",
       y="Number of Samples") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save the sample plot
ggsave(paste0(output_path,"/",sample_prefix, "_sample_missingness_histogram.png"), plot = sample_plot, width = 8, height = 6)

# Plotting Site Missingness
site_plot <- ggplot(site_missingness_data, aes(x=MissingnessThreshold, y=SiteCount)) +
  geom_bar(stat="identity", fill="salmon") +
  geom_text(aes(label=SiteCount), vjust=-0.5, color="black") +
  theme_minimal() +
  labs(title="Histogram of Site Missingness",
       x="Missingness Threshold",
       y="Number of Sites") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save the site plot
ggsave(paste0(output_path,"/",sample_prefix, "_site_missingness_histogram.png"), plot = site_plot, width = 8, height = 6)