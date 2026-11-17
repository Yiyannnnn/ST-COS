#' Map layout coordinates onto image pixels
#'
#' Scales a layout given in a y-up frame uniformly, preserving its aspect
#' ratio, centers it on the image, and returns image-pixel coordinates with
#' y increasing downward.
#'
#' @param coordinates Data frame with `x_coord` and `y_coord`.
#' @param width,height Image size in pixels.
#' @return Data frame of image-pixel coordinates.
#' @keywords internal
map_coordinates_to_image <- function(coordinates, width, height) {
  x <- coordinates$x_coord
  y <- coordinates$y_coord
  span_x <- diff(range(x, finite = TRUE))
  span_y <- diff(range(y, finite = TRUE))
  scales <- c(if (span_x > 0) width / span_x, if (span_y > 0) height / span_y)
  scale <- if (length(scales)) min(scales) else 0
  coordinates$x_coord <- (x - min(x)) * scale + (width - span_x * scale) / 2
  coordinates$y_coord <- height -
    ((y - min(y)) * scale + (height - span_y * scale) / 2)
  coordinates
}

#' Convert image-pixel coordinates to a y-up simulation frame
#'
#' @param coordinates Data frame with `x_coord` and `y_coord` in image pixels.
#' @param height Image height in pixels.
#' @return Data frame with `y_coord` measured upward from the image bottom.
#' @keywords internal
image_to_simulation_coordinates <- function(coordinates, height) {
  coordinates$y_coord <- height - coordinates$y_coord
  coordinates
}

#' Build the local ST-COS image-guided Shiny application
#'
#' Creates a local coordinate designer: simulation coordinates come from an
#' uploaded coordinate table (for example, centroids from an external
#' segmentation), image clicks, or a generated grid, hexagonal, or random
#' layout. Generated layouts are placed on an uploaded histology image without
#' changing their aspect ratio. With an image loaded, uploaded and clicked
#' coordinates are read as image pixels, and the downloaded table uses a y-up
#' frame in pixel units. Image segmentation, gene patterns, cell-type
#' composition, and dependence settings are available in the hosted
#' application ([launch_stcos_app()]).
#'
#' @return A `shiny.appobj` created by [shiny::shinyApp()].
#' @export
stcos_shiny_app <- function() {
  required <- c("shiny", "png", "jpeg")
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop(
      "Install the following packages to run the local app: ",
      paste(missing, collapse = ", "),
      "."
    )
  }

  read_image <- function(file) {
    extension <- tolower(tools::file_ext(file$name))
    image <- if (extension == "png") {
      png::readPNG(file$datapath)
    } else if (extension %in% c("jpg", "jpeg")) {
      jpeg::readJPEG(file$datapath)
    } else {
      stop("The image must be a PNG or JPEG file.")
    }
    list(image = image, width = dim(image)[2], height = dim(image)[1])
  }

  parse_coordinates <- function(path) {
    input <- utils::read.csv(path, check.names = FALSE)
    lower_names <- tolower(names(input))
    find_name <- function(candidates) {
      index <- match(candidates, lower_names, nomatch = 0L)
      index <- index[index > 0L]
      if (length(index)) names(input)[index[1L]] else NULL
    }
    x_name <- find_name(c("x_coord", "x", "col", "column", "centroid_x"))
    y_name <- find_name(c("y_coord", "y", "row", "centroid_y"))
    if (is.null(x_name) || is.null(y_name)) {
      stop("The coordinate file must contain x/y or x_coord/y_coord columns.")
    }
    coordinates <- data.frame(
      x_coord = as.numeric(input[[x_name]]),
      y_coord = as.numeric(input[[y_name]])
    )
    coordinates <- coordinates[stats::complete.cases(coordinates), , drop = FALSE]
    rownames(coordinates) <- NULL
    coordinates
  }

  ui <- shiny::fluidPage(
    shiny::tags$head(
      shiny::tags$style(shiny::HTML(
        paste(
          ".stcos-panel { border-left: 4px solid #2d7896; padding-left: 16px; }",
          ".stcos-actions .btn { margin-right: 6px; margin-bottom: 8px; }",
          ".stcos-count { font-weight: 600; color: #18324a; }"
        )
      ))
    ),
    shiny::titlePanel("ST-COS image-guided coordinate designer"),
    shiny::sidebarLayout(
      shiny::sidebarPanel(
        class = "stcos-panel",
        shiny::fileInput(
          "image_file", "Histology image",
          accept = c("image/png", "image/jpeg")
        ),
        shiny::fileInput(
          "coordinate_file", "Coordinate CSV",
          accept = c(".csv", "text/csv")
        ),
        shiny::selectInput(
          "pattern", "Generated layout",
          choices = c("Hexagonal" = "hex", "Grid" = "grid", "Random" = "random")
        ),
        shiny::fluidRow(
          shiny::column(
            6,
            shiny::numericInput("x_len", "X extent", value = 10, min = 1)
          ),
          shiny::column(
            6,
            shiny::numericInput("y_len", "Y extent", value = 10, min = 1)
          )
        ),
        shiny::fluidRow(
          shiny::column(
            6,
            shiny::numericInput("spacing", "Spot spacing", value = 1, min = 0.01)
          ),
          shiny::column(
            6,
            shiny::numericInput("n_spots", "Random spots", value = 100, min = 1)
          )
        ),
        shiny::numericInput("seed", "Random seed", value = 123, min = 1),
        shiny::tags$div(
          class = "stcos-actions",
          shiny::actionButton("generate", "Generate", icon = shiny::icon("wand-magic-sparkles")),
          shiny::actionButton("undo", "Undo", icon = shiny::icon("rotate-left")),
          shiny::actionButton("clear", "Clear", icon = shiny::icon("trash"))
        ),
        shiny::downloadButton("download_coordinates", "Download coordinates"),
        shiny::helpText(
          "With an image loaded, uploaded and clicked coordinates are image",
          "pixels; downloads use a y-up frame in pixel units."
        )
      ),
      shiny::mainPanel(
        shiny::tabsetPanel(
          shiny::tabPanel(
            "Image and coordinates",
            shiny::plotOutput("overlay", height = "650px", click = "image_click"),
            shiny::textOutput("coordinate_count", container = shiny::span)
          ),
          shiny::tabPanel(
            "Coordinate table",
            shiny::tableOutput("coordinate_table")
          )
        )
      )
    )
  )

  server <- function(input, output, session) {
    coordinates <- shiny::reactiveVal(data.frame(
      x_coord = numeric(),
      y_coord = numeric()
    ))

    image_info <- shiny::reactive({
      file <- input$image_file
      if (is.null(file)) return(NULL)
      read_image(file)
    })

    shiny::observeEvent(input$coordinate_file, {
      file <- input$coordinate_file
      shiny::req(file)
      parsed <- tryCatch(
        parse_coordinates(file$datapath),
        error = function(error) {
          shiny::showNotification(conditionMessage(error), type = "error")
          NULL
        }
      )
      if (!is.null(parsed)) coordinates(parsed)
    })

    shiny::observeEvent(input$generate, {
      generated <- generate_coordinates(
        x_len = input$x_len,
        y_len = input$y_len,
        pattern = input$pattern,
        spot_distance = if (input$pattern == "random") NULL else input$spacing,
        n_spot = if (input$pattern == "random") input$n_spots else NULL,
        seed = input$seed
      )
      generated <- as.data.frame(generated[, c("x_coord", "y_coord")])
      image <- image_info()
      if (!is.null(image)) {
        generated <- map_coordinates_to_image(generated, image$width, image$height)
      }
      coordinates(generated)
    })

    shiny::observeEvent(input$image_click, {
      click <- input$image_click
      shiny::req(click)
      current <- coordinates()
      coordinates(rbind(
        current,
        data.frame(x_coord = click$x, y_coord = click$y)
      ))
    })

    shiny::observeEvent(input$undo, {
      current <- coordinates()
      if (nrow(current)) coordinates(current[-nrow(current), , drop = FALSE])
    })

    shiny::observeEvent(input$clear, {
      coordinates(data.frame(x_coord = numeric(), y_coord = numeric()))
    })

    output$overlay <- shiny::renderPlot({
      image <- image_info()
      points <- coordinates()
      if (!is.null(image)) {
        graphics::plot.new()
        graphics::plot.window(
          xlim = c(0, image$width),
          ylim = c(image$height, 0),
          asp = 1
        )
        graphics::rasterImage(
          grDevices::as.raster(image$image),
          0, image$height, image$width, 0
        )
        graphics::box(col = "#94A3B8")
      } else if (nrow(points)) {
        graphics::plot(
          points$x_coord, points$y_coord,
          type = "n", asp = 1,
          xlab = "x", ylab = "y"
        )
      } else {
        graphics::plot.new()
        graphics::text(0.5, 0.5, "Upload an image or generate coordinates")
      }
      if (nrow(points)) {
        graphics::points(
          points$x_coord, points$y_coord,
          pch = 21, bg = "#F43F5E", col = "white", cex = 1.1
        )
      }
    }, res = 120)

    output$coordinate_count <- shiny::renderText({
      paste(format(nrow(coordinates()), big.mark = ","), "coordinates")
    })

    output$coordinate_table <- shiny::renderTable({
      utils::head(coordinates(), 100L)
    }, digits = 3, striped = TRUE, bordered = TRUE)

    output$download_coordinates <- shiny::downloadHandler(
      filename = function() "stcos_image_coordinates.csv",
      content = function(file) {
        out <- coordinates()
        image <- image_info()
        if (!is.null(image)) out <- image_to_simulation_coordinates(out, image$height)
        utils::write.csv(out, file, row.names = FALSE)
      }
    )
  }

  shiny::shinyApp(ui = ui, server = server)
}

#' Open the hosted ST-COS image-guided application
#'
#' Opens the full hosted ST-COS Shiny application, which provides interactive
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
