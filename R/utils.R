#' Start FuSViz web app
#'
#' @description The app is started with system default browser
#'
#' @param data Data to be loaded with app; the default value as NULL
#'
#' @export
FuSViz_app <- function(data = NULL) {
	DIR = system.file("app", package = "FuSViz")
	if (DIR == "") { stop("Could not find app directory. Try re-installing `FuSViz`.", call. = FALSE) }
	source(file.path(DIR, "ui.R"), local = TRUE, chdir = TRUE)
	source(file.path(DIR, "server.R"), local = TRUE, chdir = TRUE)
	shiny_app = shiny::shinyApp(ui = ui, server = server)
	shiny::runApp(shiny_app, launch.browser = TRUE, display.mode = "normal")
}

