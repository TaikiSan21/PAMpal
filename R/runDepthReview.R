#' @title Run Echo Depth Review App
#'
#' @description Runs a Shiny app to review the slant delay and esimated depths
#'   of an \linkS4class{AcousticStudy} object that has been processed with
#'   \link{calculateEchoDepth}. App allows users to select detections that
#'   should not be included in future analysis and marks them with the tag
#'   \code{keepClick=FALSE}, similar to \link{filterEchoDepths}.
#'
#' @param x an \linkS4class{AcousticStudy} object that has been processed with
#'   \link{calculateEchoDepth}
#'
#' @author Taiki Sakai \email{taiki.sakai@@noaa.gov}
#'
#' @return the object as \code{x}, with updated \code{keepClick} column
#'
#' @examples
#' # example not run because \link{calculateEchoDepth} must be run first,
#' # and it requires a large amount of data not stored in the package
#' \dontrun{
#' study <- calculateEchoDepth(study, wav='path/to/wavFiles')
#' study <- runDepthReview(x)
#' }
#'
#' @importFrom shiny fluidPage selectInput HTML fluidRow plotOutput column
#' @importFrom shiny radioButtons actionButton sliderInput div textInput
#' @importFrom shiny reactiveValues updateSelectInput observeEvent tags
#' @importFrom shiny updateSliderInput observe renderPlot showNotification
#' @importFrom shiny brushedPoints runApp shinyApp isolate onSessionEnded
#' @importFrom ggplot2 geom_point ggtitle xlab ylab theme element_text
#' @importFrom ggplot2 element_blank element_line aes scale_color_manual ylim
#' @importFrom signal butter
#'
#' @export
#'
runDepthReview <- function(x) {
    inStudy <- is.AcousticStudy(x)
    inDf <- is.data.frame(x)
    if(!inStudy && !inDf) {
        stop('Input must be AcousticStudy or dataframe')
    }
    # using SHINYDATA as a temp var name to modify with <<- later from shiny
    if(inStudy) {
        SHINYDATA <- getClickData(x)
    }
    if(inDf) {
        SHINYDATA <- x
    }
    if(!'keepClick' %in% colnames(SHINYDATA)) {
        SHINYDATA$keepClick <- !is.na(SHINYDATA$maxDepth)
    }
    SHINYDATA$maxDepth <- -1*SHINYDATA$maxDepth
    # using on.exit to return the modified data in case shiny app crashes
    # or something like that
    on.exit({
        SHINYDATA$maxDepth <- -1*SHINYDATA$maxDepth
        if(inStudy) {
            SHINYDATA <- detDataToStudy(x, SHINYDATA)
        }
        return(invisible(SHINYDATA))
    })
    #### UI ####
    ui <- fluidPage(
        selectInput('evSelect', label='Event', choices=''),
        tags$head(tags$style(HTML(".selectize-input {width: 500px;}"))),
        # brush argument will enable the brush, sends the data point information to the server side
        plotOutput(outputId = "scatterplot", brush = "plot_brush_"), # brush ID is plot_brush, brush argument enables the brush
        fluidRow(
            column(width=3,
                   radioButtons('paintFlag',
                                label='Select paint action',
                                choices=list('Remove selections'=FALSE, 'Keep selections'=TRUE)
                   ),
                   actionButton('allFalse', 'Remove all detections'),
                   radioButtons('plotValue',
                                label='Plot time delay or depth',
                                choices=list('Time delay'='maxTime', 'Depth'='maxDepth'), inline=TRUE),
                   sliderInput('yLims',
                               label='Depth limit',
                               min=-4000, max=0, value=-4000,
                               width='100%'),
                   fluidRow(div('Echogram filter values (Hz)',
                                style='font-weight:bold; font-size:14px; text-align:center'),
                            column(width=6,textInput('freqLow', label='Lower bound',
                                                     value='2000', width='100%')),
                            column(width=6, textInput('freqHigh', 'Higher bound', value='16000', width='100%'))
                   ),
                   actionButton('loadEcho',
                                label='Load echogram data (slow)'),
                   actionButton("save","Save") # when clicked saves brushed dataframe to file
            ),
            column(width=9,
                   plotOutput(outputId = 'echogram')
            )
        )
    )

    # Server code begins here ####
    server <- function(input, output, session) {
        # making the dataset reactiveValues so that any changes in mt$data later could be reflected throughout
        smf <- reactiveValues(data=SHINYDATA,
                              echoData=list())
        # this populates dropdown with event names
        updateSelectInput(inputId='evSelect', choices=unique(SHINYDATA$eventId))
        yLabels <- list(
            'maxDepth' = 'Estimated depth (m)',
            'maxTime' = 'Measured time delay (s)'
        )
        observeEvent({
            input$evSelect
            input$plotValue
        }, {
            if(input$evSelect != '') {
                yLims <- switch(
                    input$plotValue,
                    'maxDepth' = {
                        yRange <- range(smf$data$maxDepth[smf$data$eventId == input$evSelect], na.rm=TRUE)
                        yRange[1] <- round(yRange[1] -50, 0)
                        updateSliderInput(session, 'yLims',
                                          label='Depth limit', min=yRange[1], max=0, value=yRange[1],
                                          step=50)
                    },
                    'maxTime' = {
                        yRange <- range(smf$data$maxTime[smf$data$eventId == input$evSelect], na.rm=TRUE)
                        yRange[2] <- round(yRange[2] + .001, 3)
                        updateSliderInput(session, 'yLims',
                                          label='Time limit', min=0, max=yRange[2], value=yRange[2],
                                          step=.001)
                    }
                )
            }
        })
        #### scatterplot ####
        output$scatterplot <- renderPlot({
            plotData <- smf$data %>%
                dplyr::filter(.data$eventId == input$evSelect)
            g <- ggplot(plotData, aes(x = .data[['UTC']], y = .data[[input$plotValue]], col=.data[['keepClick']])) +
                geom_point(na.rm=TRUE) +
                ggtitle(input$evSelect) +
                xlab("Click time (HH:MM)") + ylab(yLabels[[input$plotValue]]) +
                theme(axis.text = element_text(size = 14), # format axis font
                      axis.title = element_text(size = 16, face = "bold"),
                      text = element_text(size=14),
                      panel.background = element_blank(),
                      axis.line = element_line(colour = "black"),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank())  +
                scale_color_manual(values=c('#F8766D', '#00BFC4'), breaks=c(FALSE, TRUE))

            if(input$plotValue == 'maxDepth') {
                g <- g + ylim(input$yLims, 0)
            }
            if(input$plotValue == 'maxTime') {
                g <- g + ylim(0, input$yLims)
            }
            g
        })

        #### load echo data ####
        observeEvent(input$loadEcho, {
            # check if we've already loaded this events data

            low <- suppressWarnings(as.numeric(input$freqLow))
            high <- suppressWarnings(as.numeric(input$freqHigh))

            # if(input$evSelect %in% names(smf$echoData)) {
            #     thisEcho <- smf$echoData[[input$evSelect]]
            # } else {
            # if we havent then load soundfiles
            showNotification('Loading wav clips...')
            thisSr <- NULL
            thisEcho <- lapply(smf$data$wavFile[smf$data$eventId == input$evSelect], function(x) {
                if(is.na(x) || !file.exists(x)) {
                    return(NULL)
                }
                wav <- readWave(x, toWaveMC = TRUE, from=0, to=.03, units='seconds')
                thisSr <<- wav@samp.rate
                wav <- wav@.Data[, 1] / 2^(wav@bit-1)
                if(!is.na(low)) {
                    lowFilt <- butter(8, 2 * low / thisSr, type='high')
                    wav <- signal::filter(lowFilt, wav)
                }
                if(!is.na(high)) {
                    high <- min(high, thisSr/2)
                    highFilt <- butter(8, 2 * high / thisSr, type='low')
                    wav <- signal::filter(highFilt, wav)
                }
                if(inherits(wav, 'ts')) {
                    wav <- unclass(wav)
                }
                wav
            })
            thisEcho <- do.call(cbind, thisEcho)
            thisEcho <- 10*log10(abs(thisEcho) + 1e-10)
            lims <- quantile(thisEcho, c(.01, .999))
            thisEcho[thisEcho < lims[1]] <- lims[1]
            thisEcho[thisEcho > lims[2]] <- lims[2]
            # store for easy access in future
            smf$echoData[[input$evSelect]] <- list(echo=thisEcho,
                                                   sr=thisSr)
            showNotification('Done loading!')
            # }
        })
        observeEvent(input$allFalse, {
            smf$data$keepClick[smf$data$eventId == input$evSelect] <- FALSE
        })
        # output from viridisLite::viridis(32)
        viridis32 <- c("#440154FF", "#470D60FF", "#48196BFF", "#482475FF", "#472E7CFF",
                       "#453882FF", "#424186FF", "#3E4B8AFF", "#3A548CFF", "#365D8DFF",
                       "#32658EFF", "#2E6D8EFF", "#2B758EFF", "#287D8EFF", "#25848EFF",
                       "#228C8DFF", "#1F948CFF", "#1E9C89FF", "#20A386FF", "#25AB82FF",
                       "#2EB37CFF", "#3ABA76FF", "#48C16EFF", "#58C765FF", "#6ACD5BFF",
                       "#7ED34FFF", "#93D741FF", "#A8DB34FF", "#BEDF26FF", "#D4E21AFF",
                       "#E9E51AFF", "#FDE725FF")
        #### echogram plot ####
        output$echogram <- renderPlot({
            if(input$evSelect %in% names(smf$echoData)) {
                toPlot <- smf$data$keepClick[smf$data$eventId == input$evSelect & !is.na(smf$data$wavFile)]
                if(any(toPlot)) {
                    echoData <- smf$echoData[[input$evSelect]]
                    echoWavs <- echoData$echo[, toPlot]
                    par(mar=c(3.1, 4.1, 3.1, 2.1))
                    image(z=t(echoWavs), x=1:ncol(echoWavs), y=(1:nrow(echoWavs))/echoData$sr*1e3,
                          xlab='', ylab='Time (ms)', main='Echogram',
                          col=viridis32, useRaster=TRUE)
                }
            }
        })

        #### observe brush ####
        observeEvent({
            input$plot_brush_ # only update on either brush or flag change
            input$paintFlag
        }, {
            df = brushedPoints(smf$data, brush = input$plot_brush_, xvar='UTC', yvar=input$plotValue, allRows = TRUE) # get column with false and true with AllRows = T
            # removing datapoints selected by brush, select values are stored as chars
            smf$data$keepClick[df$selected_ & df$eventId == input$evSelect] <- input$paintFlag == 'TRUE'
        })
        # reset brush box when we change events
        observeEvent(input$evSelect, {
            session$resetBrush('plot_brush_')
        })
        #### save progress ####
        # this is for forcing saves in case it gets buggy
        observeEvent(input$save, {
            # <<- puts it back to environment of original function call
            SHINYDATA <<- smf$data
            showNotification('Saved!')
        })
        # usually it should auto save
        onSessionEnded(function() {
            SHINYDATA <<- isolate(smf$data)
        })
    }

    # Create a Shiny app object
    app <- shinyApp(ui = ui, server = server)

    #open shiny app
    runApp(app)
}