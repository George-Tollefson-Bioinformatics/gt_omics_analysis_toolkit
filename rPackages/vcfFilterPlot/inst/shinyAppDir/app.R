library(shiny)
library(gdsfmt)
library(SeqArray)
library(ggplot2)
library(dplyr)
# default maximum file size is ~30Mb - can make this user adjustable later if running on remote server with more memory. Want to prevent crashing local computer.
options(shiny.maxRequestSize=30*1024^2)

shinyApp(ui = ui, server = server)