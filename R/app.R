#' Launch the ST-COS image-guided application
#'
#' Opens the hosted ST-COS Shiny application, which provides interactive
#' coordinate design, histology-image upload and segmentation, gene-pattern
#' configuration, cell-type composition, and gene-dependence controls.
#'
#' @param browser Logical scalar. If `TRUE`, open the application in the
#'   default web browser. If `FALSE`, return the application URL without
#'   opening a browser. The default is `interactive()`.
#'
#' @return The application URL, invisibly.
#' @export
launch_stcos_app <- function(browser = interactive()) {
  if (!is.logical(browser) || length(browser) != 1L || is.na(browser)) {
    stop("browser must be TRUE or FALSE.")
  }

  app_url <- "https://yiyanz.shinyapps.io/st-cos/"
  if (browser) utils::browseURL(app_url)
  invisible(app_url)
}
