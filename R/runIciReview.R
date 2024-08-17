#' @title Run Echo Depth Review App
#'
#' @description Runs a Shiny app 
#'
#' @param x an \linkS4class{AcousticStudy} object that has been processed with
#'   \link{calculateEchoDepth}
#' @param maxIci maximum ICI value (seconds) to display on plots
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
#' study <- runIciReview(x)
#' }
#'
#' @importFrom shiny fluidPage selectInput HTML fluidRow plotOutput column
#' @importFrom shiny actionButton stopApp icon
#' @importFrom shiny reactiveValues updateSelectInput observeEvent tags
#' @importFrom shiny renderPlot showNotification
#' @importFrom shiny  runApp shinyApp isolate onSessionEnded
#' @importFrom graphics text
#'
#' @export
#'
runIciReview <- function(x, maxIci=2.5) {
    if(!is.AcousticStudy(x)) {
        stop('"x" must be an AcousticStudy object.')
    }
    allMeasures <- getMeasures(x)
    if(!'All_ici' %in% names(allMeasures)) {
        warning('ICI data not found, calculating ICI')
        x <- calculateICI(x, 'peakTime')
    }
    if(!'ipiMax' %in% names(allMeasures)) {
        warning('IPI estimates not found')
        hasIpi <- FALSE
        avgWaves <- NULL
    } else {
        hasIpi <- TRUE
        avgWaves <- lapply(events(x), function(ev) {
            ancillary(ev)$avgWave
        })
    }
    allIci <- getClickData(x, measures=FALSE)[c('UTC', 'UID', 'Channel',
                                                'eventId', 'species', 'angle')]
    allIci <- left_join(allIci, getICI(x, 'data')[c('UID', 'Channel', 'eventId', 'ici')],
                        by=c('eventId', 'UID', 'Channel'))
    # Modify x with any updated measures 
    on.exit({
        x <- addMeasures(x, allMeasures, replace=TRUE)
        return(x)
    })
    #### UI ####
    ui <- fluidPage(
        # keyboard input ####
        tags$script('$(document).on("keydown",
                 function (e) {
                 if(e.which == 37) {
                    Shiny.onInputChange("left", new Date());
                 }
                 if(e.which == 39) {
                    Shiny.onInputChange("right", new Date());
                 }
                 });
                '),
        # top bar ####
        fluidRow(column(9, selectInput('evSelect', label='Event', choices=list('Loading...'=1))),
                 column(2, 
                        actionButton('left', label='', icon=icon('chevron-left', lib='glyphicon')),
                        actionButton('right', label='', icon=icon('chevron-right', lib='glyphicon'))),
                 column(1,
                        actionButton('stopApp', label='Stop App'))),
        tags$head(tags$style(HTML(".selectize-input {width: 500px;}"))),
        tags$head(tags$style(HTML('#left {position: absolute; top: 25px; left: 0%}'))),
        tags$head(tags$style(HTML('#right {position: absolute; top: 25px; left: 50px}'))),
        tags$head(tags$style(HTML('#stopApp {position: absolute; top: 25px;}'))),
        
        # first plot ####
        # brush argument will enable the brush, sends the data point information to the server side
        fluidRow(column(6, 
                        plotOutput(outputId = 'iciTime', click='icitime_click')),
                 column(6,
                        plotOutput(outputId = 'iciHist', click = 'hist_click'))
        ),
        # button row ####
        fluidRow(column(5, 
                        icon('ok-circle', lib='glyphicon'),
                        actionButton('saveIci', label='Save User ICI', icon=icon('chevron-up', lib='glyphicon')),
                        actionButton('saveIpi', label='Save User IPI', icon=icon('chevron-down', lib='glyphicon'))),
                 column(2,
                        actionButton('reset', label='Reset')),
                 column(5, 
                        icon('remove-circle', lib='glyphicon'),
                        actionButton('noIci', label='Remove ICI',  icon=icon('chevron-up', lib='glyphicon')),
                        actionButton('noIpi', label='Remove IPI',  icon=icon('chevron-down', lib='glyphicon')))),
        tags$head(tags$style(HTML('.glyphicon.glyphicon-remove-circle {font-size: 20px;}'))),
        tags$head(tags$style(HTML('.glyphicon.glyphicon-ok-circle {font-size: 20px}'))),
        # second plot ####
        plotOutput(outputId = "avgWave", click = "avg_click")
    )
    #### SERVER ####
    server <- function(input, output, session) {
        evList <- as.list(seq_along(events(x)))
        names(evList) <- names(events(x))
        reVals <- reactiveValues(measures=allMeasures,
                                 ici=allIci,
                                 userIpi=NA,
                                 userIci=NA,
                                 evIndex=1,
                                 thisEvent=names(evList)[1])
        
        updateSelectInput(inputId='evSelect', choices=evList)
        # event changing ####
        observeEvent(input$left, {
            if(reVals$evIndex > 1) {
                reVals$evIndex <- reVals$evIndex - 1
                updateSelectInput(inputId='evSelect', selected=as.character(reVals$evIndex))
            } else{
                showNotification('Already at first event!')
            }
        })
        observeEvent(input$right, {
            if(reVals$evIndex < length(evList)) {
                reVals$evIndex <- reVals$evIndex + 1
                updateSelectInput(inputId='evSelect', selected=as.character(reVals$evIndex))
            } else {
                showNotification('Already at last event!')
            }
        })
        observeEvent(input$evSelect, {
            reVals$evIndex <- as.numeric(input$evSelect)
            reVals$userIpi <- NA
            reVals$userIci <- NA
        })
        observeEvent(reVals$evIndex,
                     reVals$thisEvent <- names(evList)[reVals$evIndex]
        )
        # plots-ici time ####
        output$iciTime <- renderPlot({
            thisIciData <- reVals$ici[reVals$ici$eventId == reVals$thisEvent, ]
            thisIciData <- thisIciData[thisIciData$ici < maxIci, ]
            thisIci <- reVals$measures$All_ici[reVals$measures$eventId == reVals$thisEvent]
            userIci <- round(reVals$userIci, 2)
            
            par(mar=c(4,4,2,1))
            plot(x=thisIciData$UTC, y=thisIciData$ici,
                 ylim=c(0, maxIci),
                 ylab='ICI (s)',
                 xlab='UTC Time')
            lines(x=range(thisIciData$UTC), y=rep(thisIci, 2), col='red', lty=2)
            lines(x=range(thisIciData$UTC), y=rep(userIci, 2), col='blue', lty=2)
        })
        # plots-ici hist ####
        output$iciHist <- renderPlot({
            thisIciData <- reVals$ici[reVals$ici$eventId == reVals$thisEvent, ]
            thisIciData <- thisIciData[thisIciData$ici < maxIci, ]
            thisIci <- reVals$measures$All_ici[reVals$measures$eventId == reVals$thisEvent]
            userIci <- round(reVals$userIci, 2)
            par(mar=c(4,2,2,1))
            thisHist <- hist(thisIciData$ici,
                             breaks=seq(from=0, to=maxIci, by=.01),
                             main=paste0('ICI Histogram\n Modal ICI: ',
                                         round(thisIci, 2),
                                         ' User ICI: ', userIci),
                             xlab='ICI (s)',
                             ylab=''
            )
            lines(x=rep(thisIci, 2), y=c(0, 10e3), col='red', lty=2)
            lines(x=rep(userIci, 2), y=c(0, 10e3), col='blue', lty=2)
        })
        # plots-avg wave ####
        output$avgWave <- renderPlot({
            thisAvg <- avgWaves[[reVals$thisEvent]]
            if(is.null(thisAvg)) {
                plot(x=1, y=1, type='n')
                text(x=1, y=1, label='No Average Waveform')
                return()
            }
            thisMax <- which.max(thisAvg)
            origIpi <- reVals$measures$ipiMax[reVals$measures$eventId == reVals$thisEvent]
            userIpi <- round(reVals$userIpi * 1e3, 2)
            # par(mfrow=c(1,2))
            plot(x=(seq_along(thisAvg)-thisMax)/attr(thisAvg, 'sr') * 1e3,
                 y=thisAvg,
                 type='l',
                 main=paste0('Average Waveform\n Orig IPI: ', round(origIpi*1e3, 2),
                             ' User IPI: ', userIpi),
                 xlab='IPI Est (ms)',
                 ylab='')
            lines(x=rep(origIpi * 1e3, 2), y=c(-1, 1), col='red', lty=2)
            lines(x=rep(userIpi, 2), y=c(-1, 1), col='blue', lty=2)
            
        })
        # clicking on plots ####
        observeEvent(input$hist_click, {
            reVals$userIci <- input$hist_click$x
        })
        observeEvent(input$icitime_click, {
            reVals$userIci <- input$icitime_click$y
        })
        observeEvent(input$avg_click, {
            reVals$userIpi <- input$avg_click$x / 1e3
            
        })
        onSessionEnded(function() {
            allMeasures <<- isolate(reVals$measures)
        })
        observeEvent(input$stopApp, {
            allMeasures <<- isolate(reVals$measures)
            stopApp()
        })
        # save/na buttons ####
        observeEvent(input$saveIci, {
            reVals$measures$All_ici[reVals$measures$eventId == reVals$thisEvent] <- reVals$userIci
        })
        observeEvent(input$saveIpi, {
            reVals$measures$ipiMax[reVals$measures$eventId == reVals$thisEvent] <- reVals$userIpi
        })
        observeEvent(input$noIci, {
            reVals$measures$All_ici[reVals$measures$eventId == reVals$thisEvent] <- NA
        })
        observeEvent(input$noIpi, {
            reVals$measures$ipiMax[reVals$measures$eventId == reVals$thisEvent] <- NA
        })
        observeEvent(input$reset, {
            thisAvg <- avgWaves[[reVals$thisEvent]]
            avgEcho <- findEchoTimes(wav=thisAvg, sr=attr(thisAvg, 'sr'), 
                                     clipLen=length(thisAvg)/attr(thisAvg, 'sr'), plot=FALSE)
            thisIciData <- reVals$ici[reVals$ici$eventId == reVals$thisEvent, ]
            thisIci <- calcIciMode(thisIciData$ici, c(0, maxIci))
            reVals$measures$All_ici[reVals$measures$eventId == reVals$thisEvent] <- thisIci
            reVals$measures$ipiMax[reVals$measures$eventId == reVals$thisEvent] <- avgEcho$time[1]
        })
    }
    app <- shinyApp(ui, server)
    runApp(app)
}

# ipiData <- readRDS('../Data/PmICI/HB2102_IpiTestStudy.rds')
### CURRENT PROBLEMS
## Clicking doesnt work with par mfrow, interacts with second (last) plot
## Wont be able to interact on both, will have to not do the par way.
## can google if theres a way to work click input across multi plot
## Right now the saving part doesnt actually do anything or save it
# new <- runIciReview(ipiData)
