#' @title Plot Spectrogram or Cepstrogram
#'
#' @description Plots either a spectrogram or cepstrogram and also overlays
#'   whistle or cepstral contours from the binary files
#'
#' @param x an \linkS4class{AcousticStudy} object
#' @param evNum if \code{x} is a study, the event index number to calculate the average
#'   spectra for. Note that this is the index in the order that they appear in the
#'   \linkS4class{AcousticStudy} object, not the actual event number. Alternatively
#'   full event names can be used
#' @param start start time of the plot, either POSIXct or seconds from
#'   the start of the recording
#' @param end end time of the plot, either POSIXct or seconds from the start
#'   of the recording.
#' @param channel channel to plot
#' @param wl window length of FFT to use for creating plot
#' @param hop hop value of FFT to use for creating plot, either as a percentage
#'   of \code{wl} or number of samples
#' @param mode one of \code{'spec'} or \code{'ceps'} to plot either spectrogram
#'   or cepstrogram
#' @param detections vector containing any of \code{'whistle'}, \code{'click'},
#'   and/or \code{'cepstrum'} indicating which detections to overlay on the plot
#' @param sr sample rate
#'
#' @return nothing, just plots
#'
#' @author Taiki Sakai \email{taiki.sakai@@noaa.gov}
#'
#' @examples
#' 
#' data(exStudy)
#' recs <- system.file('extdata', 'Recordings', package='PAMpal')
#' exStudy <- addRecordings(exStudy, folder=recs, log=FALSE, progress=FALSE)
#' # No detections will appear on plot because included recordings are heavily decimated
#' plotGram(exStudy)
#' 
#' @importFrom dplyr bind_rows
#' @importFrom tuneR readWave downsample
#' @importFrom grDevices gray.colors
#' @importFrom graphics points
#'
#' @export
#'
plotGram <- function(x, evNum=1,  start=0, end=NULL, channel=1, wl=512, hop=.25, mode=c('spec', 'ceps'),
                     detections = c('whistle', 'click', 'cepstrum'), sr=NULL) {
    # Needs to be one event at a time
    x <- x[evNum]
    fileExists <- checkRecordings(x)
    thisRec <- files(x)$recordings
    thisRec <- thisRec[thisRec$db == files(x[[1]])$db, ]
    if(hop <=1) {
        hop <- hop * wl
    }
    dets <- bind_rows(lapply(getDetectorData(x), function(d) {
        list(UTC = d$UTC, UID = d$UID)
    }))

    fileCheck <- checkIn(min(dets$UTC), thisRec)
    if(is.na(fileCheck)) return(x)
    if(is.null(end)) {
        end <-thisRec$length[fileCheck]
    }
    if(inherits(start, 'POSIXct')) {
        start <- as.numeric(difftime(start, thisRec$start[fileCheck], units='secs'))
    }
    if(inherits(end, 'POSIXct')) {
        end <- as.numeric(difftime(end, thisRec$start[fileCheck], units='secs'))
    }
    wav <- readWave(thisRec$file[fileCheck], from=start, to=end, units='seconds', toWaveMC = TRUE)
    if(is.null(sr)) {
        sr <- wav@samp.rate
    }
    if(sr != wav@samp.rate) {
        wav <- downsample(wav, sr)
    }
    timeStart <- thisRec$start[fileCheck] + start
    timeEnd <- timeStart - start + end

    if(channel > ncol(wav@.Data)) {
        stop('Specified channel is not present in wav file.', call.=FALSE)
    }
    wav <- wav@.Data[, channel]
    nSlices <- ceiling((length(wav) - wl)/(hop)) + 1
    slices <- seq(from=1, by=hop, length.out=nSlices)
    data <- apply(as.matrix(slices), 1, function(s) {
        thisWav <- wav[s:(s+wl-1)]
        thisGram <- myGram(thisWav, mode=mode, wl=wl, channel=channel, sr=sr)
        thisGram[, 2]
    })

    x <- filter(x, UTC <= timeEnd, UTC >= timeStart)

    yAxis <- myGram(wav[1:wl], mode=mode, wl=wl, channel=channel, sr=sr)[, 1]
    switch(match.arg(mode),
           'spec' = {
               title <- 'Spectrogram'
               yName = 'Frequency (kHz)'
               yAxis <- yAxis / 1e3
               plotMat <- data
           },
           'ceps' = {
               title <- 'Cepstrogram'
               yName = 'ICI (ms)'
               yAxis <- yAxis * 1e3
               plotMat <- 20*log10(abs(data))
           }
    )
    xName <- paste0('Seconds after ', timeStart)
    image(t(plotMat), col = gray.colors(64, start=1, end=0), xaxt='n', yaxt='n',
          xlab = xName, ylab=yName)
    tPretty <- pretty((start:end) - start, n=5)
    # tLabs <- files(x)$recordings$start[fileCheck] + tPretty
    tLabs <- tPretty
    tLocs <- (tPretty)/(end-start)
    axis(1, at=tLocs, labels=tLabs)

    yPretty <- pretty(yAxis, n=5)
    axis(2, at=yPretty/max(yAxis), labels = yPretty)
    title(main=title)

    if('cepstrum' %in% detections) {
        plotCeps(x, timeStart, timeEnd)
    }
    if('whistle' %in% detections) {
        plotWhistle(x, timeStart, timeEnd)
    }
    if('click' %in% detections) {
        plotClick(x, timeStart, timeEnd, yMin=0, yMax = sr / 2 / 1e3)
    }
    invisible(data)
}

plotCeps <- function(x, tMin, tMax) {
    UIDs <- getCepstrumData(x)$UID

    if(is.null(UIDs) ||
       length(UIDs) == 0) {
        return(NULL)
    }
    binData <- getBinaryData(x, UIDs, type='whistle')
    settings <- getPamFft(binData)
    settings$sr <- settings$sr * 2
    settings$wl <- settings$wl * 2
    for(i in seq_along(binData)) {
        plotOneCeps(binData[[i]], hop=settings$hop, sr=settings$sr,
                    yMin=0, yMax=settings$wl/2, tMin=tMin, tMax=tMax)
    }
}
plotWhistle <- function(x, tMin, tMax) {
    UIDs <- getWhistleData(x)$UID
    if(is.null(UIDs) ||
       length(UIDs) == 0) {
        return(NULL)
    }
    binData <- getBinaryData(x, UIDs, type='whistle')
    settings <- getPamFft(binData)
    for(i in seq_along(binData)) {
        plotOneCeps(binData[[i]], hop=settings$hop, sr=settings$sr,
                    yMin=0, yMax=settings$wl/2, tMin=tMin, tMax=tMax)
    }
}

plotOneCeps <- function(x, hop, sr, yMin=0, yMax, tMin, tMax) {
    xVals <- x$date + seq(from=0, by=hop/sr, length.out=x$nSlices)
    yVals <- x$contour

    xPlot <- scaleToOne(xVals, tMin, tMax)
    yPlot <- scaleToOne(yVals, yMin, yMax)
    lines(x=xPlot, y=yPlot, col='blue', lwd=2)
}

plotClick <- function(x, tMin, tMax, yMin, yMax) {
    clickData <- getClickData(x)[, c('UTC', 'peak')]
    if(is.null(clickData) ||
       nrow(clickData) == 0) {
        return(NULL)
    }
    xVals <- clickData$UTC
    yVals <- clickData$peak
    xPlot <- scaleToOne(xVals, tMin, tMax)
    yPlot <- scaleToOne(yVals, yMin, yMax)
    points(x=xPlot, y=yPlot, pch=1, col='blue')
}
scaleToOne <- function(vals, min, max) {
    vals <- as.numeric(vals)
    min <- as.numeric(min)
    max <- as.numeric(max)
    (vals - min) / (max-min)
}
