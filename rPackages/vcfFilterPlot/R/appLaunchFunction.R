runMyApp <- function() {
    appDir <- system.file("shinyAppDir", package = "vcfFilterPlot")
    if (appDir == "") {
        stop("Shiny app not found in the package")
    }
    shiny::runApp(appDir)
}