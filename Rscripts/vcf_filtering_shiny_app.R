library(shiny)
library(gdsfmt)
library(SeqArray)
library(ggplot2)
library(dplyr)

options(shiny.maxRequestSize=30*1024^2)

#' Open a GDS file from a VCF file
#'
#' This function checks if a GDS file already exists for a given VCF file path.
#' If the GDS file does not exist, it converts the VCF file into GDS format
#' and then opens the GDS file.
#'
#' @param vcfFilePath Path to the VCF file.
#' @param gdsFilePath Path where the GDS file will be saved or is already saved.
#'
#' @return A GDS file object.
#' @export
#'
#' @examples
#' vcfFile <- "path/to/your.vcf.gz"
#' gdsFile <- "path/to/your.gds"
#' gdsObj <- openGDSFile(vcfFile, gdsFile)
openGDSFile <- function(vcfFilePath, gdsFilePath) {
  # Check if the GDS file already exists
  if (!file.exists(gdsFilePath)) {
    # If the GDS file does not exist, convert the VCF to GDS format
    message("Converting VCF to GDS format...")
    SeqArray::seqVCF2GDS(vcfFilePath, gdsFilePath, verbose = TRUE)
  }
  
  # Now open the GDS file
  message("Opening GDS file...")
  gds <- SeqArray::seqOpen(gdsFilePath)
  
  # Return the GDS object for further use
  return(gds)
}

#' Shiny UI for VCF Missingness Analysis
#'
#' This UI sets up the user interface for the VCF Missingness Analysis Shiny application.
#' It includes inputs for uploading a VCF file, setting thresholds for sample and site missingness, 
#' and for QUAL scores. It also provides action buttons to run the analysis.
#' 
#' @return A Shiny UI object.
#' @keywords internal
ui <- fluidPage(
  titlePanel("VCF Missingness Analysis"),
  sidebarLayout(
    sidebarPanel(
      fileInput("vcfFile", "Upload VCF File", accept = ".vcf.gz"),
      textInput("outputFolder", "Output Folder", value = getwd()),
      
      # Slider and text input for Sample Missingness Threshold
      fluidRow(
        column(6, sliderInput("sampleMissThresh", "Sample Missingness Threshold (%)", 
                              min = 0, max = 100, value = 50)),
        column(6, textInput("sampleMissText", "Enter Sample Threshold (%)", value = "50"))
      ),
      
      # Slider and text input for Site Missingness Threshold
      fluidRow(
        column(6, sliderInput("siteMissThresh", "Site Missingness Threshold (%)", 
                              min = 0, max = 100, value = 90)),
        column(6, textInput("siteMissText", "Enter Site Threshold (%)", value = "90"))
      ),
      
      # Fefault QUAL threshold is 100 which is  0.00000001% chance of incorrect variant call
      # QUAL score is usually given as a Phred score, which represents the probability of an error on a logarithmic scale.
      # There is no set upper limit, using 1000000 since I see range is ~0-1 million in test data and that some high site missingness variants drop off around 500k threshold
      fluidRow(
        column(6, sliderInput("qualThresh", "Variant QUAL Threshold",
                              min = 0, max = 1000000, value = 100)),
        column(6, textInput("qualText", "Enter QUAL Threshold", value = "100"))
      ),
      
      actionButton("runAnalysis", "Run Analysis")
    ),
    mainPanel(
      plotOutput("sampleMissPlot"),
      plotOutput("siteMissPlot")
    )
  )
)

#' Shiny Server Function for VCF Missingness Analysis
#'
#' Handles the server-side logic for the VCF Missingness Analysis Shiny application.
#' It processes user inputs, performs filtering based on set thresholds, 
#' and generates plots for sample and site missingness as well as QUAL score distributions.
#'
#' @param input, output, session Reactive elements from Shiny UI.
#' @keywords internal
server <- function(input, output, session) {
  
# Define sliders and text boxes, and make them update each other automatically when one is changed.
  observeEvent(input$updatePlot, {
    updateSliderInput(session, "sampleMissThresh", value = as.numeric(input$sampleMissText))
  })
  
  observeEvent(input$sampleMissThresh, {
    updateTextInput(session, "sampleMissText", value = as.character(input$sampleMissThresh))
  })
  
  observeEvent(input$sampleMissText, {
    if (is.numeric(as.numeric(input$sampleMissText))) {
      updateSliderInput(session, "sampleMissThresh", value = as.numeric(input$sampleMissText))
    }
  })
  
  observeEvent(input$qualText, {
    if (is.numeric(as.numeric(input$qualText))) {
      updateSliderInput(session, "qualThresh", value = as.numeric(input$qualText))
    }
  })
  
  observeEvent(input$qualThresh, {
    updateTextInput(session, "qualText", value = as.character(input$qualThresh))
  })
  
  # Use reactive expression to monitor changes in both slider and text input
  reactiveThreshold <- reactive({
    # Update slider based on text input
    if (!is.null(input$sampleMissText) && input$sampleMissText != "") {
      if (is.numeric(as.numeric(input$sampleMissText))) {
        threshold <- as.numeric(input$sampleMissText)
        updateSliderInput(session, "sampleMissThresh", value = threshold)
      }
    } else {
      threshold <- input$sampleMissThresh
    }
    threshold
  })
  
  observeEvent(input$runAnalysis, {
    
    # Ensure a file has been uploaded
    req(input$vcfFile)
    inFile <- input$vcfFile$datapath
    
    # Determine the output directory and ensure it exists
    outputDir <- normalizePath(input$outputFolder)
    if (!dir.exists(outputDir)) {
      dir.create(outputDir, recursive = TRUE)
    }
    
    # Construct the GDS file path based on the VCF file name
    gdsFilePath <- file.path(outputDir, paste0(tools::file_path_sans_ext(basename(inFile)), ".gds"))
    
    # Call openGDSFile to handle the GDS file opening/creation
    gds <- openGDSFile(inFile, gdsFilePath)
    
    # define variant and sample ids once
    all_variant_ids <- seqGetData(gds, "variant.id")
    all_sample_ids <- seqGetData(gds, "sample.id")
    
    # filter for QUAL:
    qualScores <- seqGetData(gds, "annotation/qual")
    qual_variants_to_include <- which(qualScores >= input$qualThresh)
    qual_variant_ids_to_include <- all_variant_ids[qual_variants_to_include]
    seqSetFilter(gds, variant.id=qual_variant_ids_to_include)
    
    # Later dev: GQ filter will need to be done on a per sample basis so will be more nuanced (apply filter to each sample across all variants - vectorize this)
    # # Filtering for genotype quality first - will not display in sample and site missingess calc tables otherwise currently
    # gqData <- seqGetData(gds, "annotation/format/GQ")
    # variants_failing_gq_to_exclude <- apply(gqData, 1, function(gqRow) any(gqRow < input$gqThreshold))
    # variant_passing_gq_to_include <- !variants_to_exclude
    # variant_ids_passing_gq_to_include <- all_variant_ids[variant_passing_gq_to_include]
    # seqSetFilter(gds, variant.id=variant_ids_passing_gq_to_include)
  
    # Filtering for site missingness
    site_missingness <- seqMissing(gds, per.variant = TRUE)
    variants_to_include <- which(site_missingness <= input$siteMissThresh / 100)
    variant_ids_to_include <- all_variant_ids[variants_to_include]
    seqSetFilter(gds, variant.id=variant_ids_to_include)
    
    # Filtering for sample missingness
    sample_missingness <- seqMissing(gds, per.variant = FALSE)
    samples_to_include <- which(sample_missingness <= input$sampleMissThresh / 100)
    sample_ids_to_include <- all_sample_ids[samples_to_include]
    seqSetFilter(gds, sample.id=sample_ids_to_include)
    
    # Set up sample missingness plotting:
    sample_missingness_df <- data.frame(Missingness = sample_missingness[samples_to_include])
    sample_bins <- seq(0, 1, by=0.1)
    sample_labels <- paste0(seq(0, 90, by=10), "%-", seq(10, 100, by=10), "%")
    sample_missingness_df$Bin <- cut(sample_missingness_df$Missingness, sample_bins, labels = sample_labels, include.lowest = TRUE)
    # Create a data frame with all possible sample_bins to ensure all sample_bins are included in the plot
    sample_all_bins_df <- data.frame(Bin = sample_labels)
    # Counting occurrences per bin
    sample_count_df <- sample_missingness_df %>% group_by(Bin) %>% summarise(n = n())
    # Merging with all possible bins to ensure zeros are included
    sample_plot_df <- merge(sample_all_bins_df, sample_count_df, by = "Bin", all.x = TRUE)
    sample_plot_df$n[is.na(sample_plot_df$n)] <- 0  # Replacing NA with 0
    
    # Set up site missingness plotting:
    site_missingness_df <- data.frame(Missingness = site_missingness[variants_to_include])
    site_bins <- seq(0, 1, by=0.1)
    site_labels <- paste0(seq(0, 90, by=10), "%-", seq(10, 100, by=10), "%")
    site_missingness_df$Bin <- cut(site_missingness_df$Missingness, site_bins, labels = site_labels, include.lowest = TRUE)
    
    # Creating a complete set of bins for site missingness
    site_all_bins_df <- data.frame(Bin = site_labels)
    
    # Counting occurrences per bin for site missingness
    site_count_df <- site_missingness_df %>% dplyr::group_by(Bin) %>% dplyr::summarise(n = n())
    
    # Merging with all possible bins to ensure zeros are included for site missingness
    site_plot_df <- merge(site_all_bins_df, site_count_df, by = "Bin", all.x = TRUE)
    site_plot_df$n[is.na(site_plot_df$n)] <- 0  # Replacing NA with 0
    
    # plot sample missingness
    output$sampleMissPlot <- renderPlot({
      
      ggplot(sample_plot_df, aes(x = Bin, y = n)) +
        geom_bar(stat = "identity", fill = "steelblue") +
        geom_vline(xintercept = as.numeric(input$sampleMissThresh) / 10 + 0.5, color = "red", linetype = "dashed", linewidth = 1) +
        labs(title = "Sample Missingness", x = "Missingness Range", y = "Count") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    })
    
    # plot site missingness
    output$siteMissPlot <- renderPlot({
      ggplot(site_plot_df, aes(x = Bin, y = n)) +
        geom_bar(stat = "identity", fill = "salmon") +
        geom_vline(xintercept = as.numeric(input$siteMissThresh) / 10 + 0.5, color = "red", linetype = "dashed", linewidth = 1) +
        labs(title = "Site Missingness", x = "Missingness Range", y = "Count") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    })
    
    # Close the GDS file after analysis
    seqClose(gds)
  })
}

shinyApp(ui = ui, server = server)