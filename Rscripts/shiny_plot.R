library(shiny)
library(ggplot2)

# Define UI
ui <- fluidPage(
  titlePanel("Sample and Site Missingness Visualization"),
  sidebarLayout(
    sidebarPanel(
      fileInput('sampleFile', 'Choose Sample Missingness File',
                accept = c('text/tsv', 'text/tab-separated-values,text/plain', '.tsv')),
      fileInput('siteFile', 'Choose Site Missingness File',
                accept = c('text/tsv', 'text/tab-separated-values,text/plain', '.tsv')),
      sliderInput('sampleSlider', 'Sample Missingness Threshold', min = 0, max = 100, value = 50),
      sliderInput('siteSlider', 'Site Missingness Threshold', min = 0, max = 100, value = 50),
      downloadButton('downloadSample', 'Download Sample Plot'),
      downloadButton('downloadSite', 'Download Site Plot')
    ),
    mainPanel(
      plotOutput("samplePlot"),
      plotOutput("sitePlot")
    )
  )
)

# Define server logic
server <- function(input, output) {
  sampleData <- reactive({
    inFile <- input$sampleFile
    if(is.null(inFile))
      return(NULL)
    read.table(inFile$datapath, header = TRUE, sep = "\t")
  })
  
  siteData <- reactive({
    inFile <- input$siteFile
    if(is.null(inFile))
      return(NULL)
    read.table(inFile$datapath, header = TRUE, sep = "\t")
  })

  output$samplePlot <- renderPlot({
    data <- sampleData()
    if(is.null(data))
      return(NULL)
    ggplot(data, aes(x=as.factor(MissingnessThreshold), y=SampleCount)) +
      geom_bar(stat="identity", fill="skyblue") +
      geom_text(aes(label=SampleCount), vjust=-0.5, color="black") +
      geom_vline(xintercept = as.numeric(as.character(input$sampleSlider))/10 + 0.5, color = "red") + # Adjust for percentage bin labels
      theme_minimal() +
      labs(title="Histogram of Sample Missingness",
           x="Missingness Threshold",
           y="Number of Samples") +
      scale_x_discrete(breaks=1:10, labels=paste0(seq(0, 90, by=10), "%-", seq(10, 100, by=10), "%")) + # Custom x-axis labels
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  })
  
  output$sitePlot <- renderPlot({
    data <- siteData()
    if(is.null(data))
      return(NULL)
    ggplot(data, aes(x=as.factor(MissingnessThreshold), y=SiteCount)) +
      geom_bar(stat="identity", fill="salmon") +
      geom_text(aes(label=SiteCount), vjust=-0.5, color="black") +
      geom_vline(xintercept = as.numeric(as.character(input$siteSlider))/10 + 0.5, color = "red") + # Adjust for percentage bin labels
      theme_minimal() +
      labs(title="Histogram of Site Missingness",
           x="Missingness Threshold",
           y="Number of Sites") +
      scale_x_discrete(breaks=1:10, labels=paste0(seq(0, 90, by=10), "%-", seq(10, 100, by=10), "%")) + # Custom x-axis labels
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  })

  # Add download handlers
  output$downloadSample <- downloadHandler(
    filename = function() { paste("sample_missingness_histogram", Sys.Date(), ".png", sep="") },
    content = function(file) {
      data <- sampleData()
      ggsave(file, plot = output$samplePlot(), width = 8, height = 6)
    }
  )
  
  output$downloadSite <- downloadHandler(
    filename = function() { paste("site_missingness_histogram", Sys.Date(), ".png", sep="") },
    content = function(file) {
      ggsave(file, plot = output$sitePlot(), width = 8, height = 6)
    }
  )
  
}

# Run the application
shinyApp(ui = ui, server = server)