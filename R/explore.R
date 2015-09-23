#' A GUI tool for exploring purity results
#' 
#' @export
explore = function() {
  appDir = system.file("shiny", "explorer", package = "supure")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `supure`.", call. = FALSE)
  }
  
  shiny::runApp(appDir, display.mode = "normal")
}