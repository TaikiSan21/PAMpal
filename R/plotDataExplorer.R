#' @title Explore Data in an Interactive Shiny Plot
#'
#' @description Runs an interactive Shiny plot of detector data. Allows user to
#'   choose which numeric data to plot, and will allow user to both color and
#'   facet the plot by event number, detector name, or species
#'
#' @param x data to plot, can be an \code{AcousticStudy}, \code{AcousticEvent},
#'   data.frame or a list of \code{AcousticEvent} objects
#'
#' @return nothing, just plots
#'
#' @author Taiki Sakai \email{taiki.sakai@@noaa.gov}
#'
#' @examples
#'
#' data(exStudy)
#'
#' if(interactive()) plotDataExplorer(exStudy)
#' if(interactive()) plotDataExplorer(exStudy)
#'
#' @importFrom shiny fluidPage plotOutput fluidRow column selectInput
#' @importFrom shiny updateSelectInput reactiveVal observeEvent renderPlot
#' @importFrom shiny shinyApp runApp
#' @importFrom ggplot2 ggplot geom_density facet_wrap
#' @export
#'
plotDataExplorer <- function(x) {
    dets <- getDetectorData(x)

    ui <- fluidPage(
        plotOutput(outputId = 'mainPlot', height='600px'),
        fluidRow(
            #controls
            column(width=6,
                   selectInput('plotDetector', label='Detector to Plot', choices=names(dets)),
                   selectInput('plotValue', label='Value to Plot', choices='')
            ),
            column(width=6,
                   selectInput('plotColor', label='Color by...',
                               choices=c('None', 'eventId', 'detectorName', 'species'),
                               selected='None'),
                   selectInput('plotFacet', label='Facet by...',
                               choices=c('None', 'eventId', 'detectorName', 'species'),
                               selected='None')
            )
        )
    )

    server <- function(input, output, session) {
        noPlotCols <- c('UID', 'Channel', 'BinaryFile', 'eventId',
                        'detectorName', 'db', 'species')
        updateSelectInput(inputId='plotDetector', choices=names(dets))
        plotData <- reactiveVal(NULL)
        colorCol <- reactiveVal(NULL)
        facetCol <- reactiveVal(NULL)
        observeEvent(input$plotDetector, {
            plotData(dets[[input$plotDetector]])
            cols <- colnames(plotData())
            cols <- cols[!cols %in% noPlotCols]
            updateSelectInput(inputId='plotValue', choices=cols)
        })
        observeEvent(input$plotColor,
                     if(input$plotColor == 'None') {
                         colorCol(NULL)
                     } else {
                         colorCol(input$plotColor)
                     }
        )
        observeEvent(input$plotFacet,
                     if(input$plotFacet == 'None') {
                         facetCol(NULL)
                     } else {
                         facetCol(input$plotFacet)
                     }
        )
        output$mainPlot <- renderPlot({
            # cat(str(plotData()))
            if(is.null(plotData()) ||
               input$plotValue == '' ||
               !input$plotValue %in% colnames(plotData())) {
                return(ggplot())
            }

            if(is.null(colorCol())) {
                plotOut <- plotData() %>%
                    ggplot(aes(x=.data[[input$plotValue]]))
            } else {
                plotOut <- plotData() %>%
                    ggplot(aes(x=.data[[input$plotValue]],
                               col=.data[[colorCol()]]))
            }
            plotOut <- plotOut +
                geom_density()
            if(!is.null(facetCol())) {
                plotOut <- plotOut +
                    facet_wrap(facets=facetCol(), scales='free')
            }
            plotOut

        })

    }
    app <- shinyApp(ui, server)

    runApp(app)
}

# plotDataExplorer <- function(x, callType=NULL, maxCategories=15) {
#     x <- getDetectorData(x)
#     if(is.list(x)) {
#         if(length(x) == 1) {
#             x <- x[[1]]
#         } else if(length(x) > 1) {
#             if(!is.null(callType) &&
#                callType %in% names(x)) {
#                 x <- x[[callType]]
#             } else {
#                 choice <- menu(choices = names(x), title = 'What call type would you like to explore?')
#                 x <- x[[choice]]
#             }
#         } else { # case length 0 data, shouldn't happen
#             stop('No data found')
#         }
#     }
#     if(!is.data.frame(x)) { # this shouldnt be possible
#         stop('Could not convert "x" to a data frame.')
#     }
#     charCols <- colnames(x)[sapply(x, function(c) is.character(c) || is.factor(c))]
#     nColors <- sapply(charCols, function(c) {
#         length(unique(x[[c]]))
#     })
#     charCols <- charCols[nColors <= maxCategories]
#     numCols <- colnames(x)[sapply(x, function(c) is.numeric(c) || inherits(c, 'POSIXct'))]
#
#     dPick <- myPicker(numCols, label = 'Data')
#     cPick <- myPicker(charCols, none=TRUE, label = 'Color')
#     facetPick <- myPicker(charCols, none=TRUE, label = 'Facet By')
#
#     manipulate({
#         ggplot(x) +
#             geom_density(aes_string(x=dataPicker, col=colPicker)) +
#             facet_wrap(facets = fCheck)
#
#     },
#     dataPicker =dPick,
#     colPicker = cPick,
#     fCheck = facetPick)
# }
#
# globalVariables(c('dataPicker', 'colPicker', 'fCheck'))
#
# myPicker <- function(x, none=FALSE, label) {
#     vals <- as.list(x)
#     if(none) {
#         vals <- c('None' = list(NULL), vals)
#         x <- c('None', x)
#     }
#     names(vals) <- x
#
#     myPick <- list(type = 1,
#                    choices = x,
#                    values = vals,
#                    initialValue = x[1],
#                    label = label)
#     class(myPick) <- 'manipulator.picker'
#     myPick
# }
