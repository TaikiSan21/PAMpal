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
#'   the start of the event
#' @param end end time of the plot, either POSIXct or seconds from the start
#'   of the event.
#' @param channel channel to plot
#' @param wl window length of FFT to use for creating plot
#' @param hop hop value of FFT to use for creating plot, either as a percentage
#'   of \code{wl} or number of samples
#' @param mode one of \code{'spec'} or \code{'ceps'} to plot either spectrogram
#'   or cepstrogram
#' @param detections vector containing any of \code{'cepstrum'}, \code{'click'},
#'   and/or \code{'whistle'} indicating which detections to overlay on the plot
#' @param detCol vector containing colors to use for plotting detections. Order matches
#'   order of detections (default alphabetical - cepstrum, click, whistle)
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
#' exStudy <- updateFiles(exStudy,
#'                        bin=system.file('extdata', 'Binaries', package='PAMpal'),
#'                        db = system.file('extdata', 'Example.sqlite3', package='PAMpal'))
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
plotGram <- function(x, evNum=1,  start=NULL, end=NULL, channel=1, wl=512, hop=.25, mode=c('spec', 'ceps'),
                     detections = c('cepstrum', 'click', 'whistle'), detCol=c('red', 'blue', 'blue'), sr=NULL) {
    # Needs to be one event at a time
    x <- x[evNum]
    thisRec <- checkRecordings(x)
    # thisRec <- files(x)$recordings
    thisRec <- thisRec[basename(thisRec$db) == basename(files(x[[1]])$db), ]
    if(hop <=1) {
        hop <- hop * wl
    }
    dets <- bind_rows(lapply(getDetectorData(x), function(d) {
        list(UTC = d$UTC, UID = d$UID)
    }))

    if(is.null(start)) {
        start <- min(dets$UTC)
    }
    if(is.null(end)) {
        end <- max(dets$UTC)
    }
    if(is.numeric(start)) {
        start <- min(dets$UTC) + start
    }
    if(is.numeric(end)) {
        end <- min(dets$UTC) + end
    }
    # startIn <- checkIn(start, thisRec)
    # if(is.na(startIn)) {
    #     stop('Could not find valid matching recording files for these detections')
    # }
    # if(inherits(start, 'POSIXct')) {
    #     start <- as.numeric(difftime(start, thisRec$start[fileCheck], units='secs'))
    # }
    # if(inherits(end, 'POSIXct')) {
    #     end <- as.numeric(difftime(end, thisRec$start[fileCheck], units='secs'))
    # }
    # wav <- readWave(thisRec$file[fileCheck], from=start, to=end, units='seconds', toWaveMC = TRUE)
    startBuff <- as.numeric(difftime(min(dets$UTC), start, units='secs'))
    endBuff <- as.numeric(difftime(end, max(dets$UTC), units='secs'))
    wav <- getClipData(x, buffer=c(startBuff, endBuff), mode='event', channel=channel, verbose=FALSE, progress=FALSE)[[1]]
    if(is.null(sr)) {
        sr <- wav@samp.rate
    }
    if(sr > wav@samp.rate) {
        stop('Chosen sampling rate is higher than the wav files sample rate, cannot upsample.')
    }
    if(sr != wav@samp.rate) {
        wav <- downsample(wav, sr)
    }
    mode <- match.arg(mode)
    # timeStart <- thisRec$start[fileCheck] + start
    # timeEnd <- timeStart - start + end
    timeStart <- start
    timeEnd <- timeStart + length(wav) / sr

    if(channel > ncol(wav@.Data)) {
        stop('Specified channel is not present in wav file.', call.=FALSE)
    }
    wav <- wav@.Data[, channel]

    data <- wavToGram(wav, wl=wl, hop=hop, sr=sr, mode=mode)

    x <- filter(x, UTC <= timeEnd, UTC >= timeStart)

    yAxis <- myGram(wav[1:wl], mode=mode, wl=wl, channel=channel, sr=sr)[, 1]
    switch(mode,
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
    image(plotMat, col = gray.colors(64, start=1, end=0), xaxt='n', yaxt='n',
          xlab = xName, ylab=yName, useRaster=TRUE)
    timeNum <- c(as.numeric(timeStart), as.numeric(timeEnd)) - as.numeric(timeStart)
    tPretty <- pretty(timeNum, n=5)
    # tLabs <- files(x)$recordings$start[fileCheck] + tPretty
    tLabs <- tPretty
    tLocs <- (tPretty)/(diff(timeNum))
    axis(1, at=tLocs, labels=tLabs)

    yPretty <- pretty(yAxis, n=5)
    axis(2, at=yPretty/max(yAxis), labels = yPretty)
    title(main=title)
    if(length(detCol) < length(detections)) {
        detCol <- rep(detCol, length.out=length(detections))
    }
    if('cepstrum' %in% detections) {
        plotCeps(x, timeStart, timeEnd, col=detCol[detections == 'cepstrum'])
    }
    if('whistle' %in% detections) {
        plotWhistle(x, timeStart, timeEnd, plotSr = sr, col=detCol[detections == 'whistle'])
    }
    if('click' %in% detections) {
        plotClick(x, timeStart, timeEnd, yMin=0, yMax = sr / 2 / 1e3, col=detCol[detections == 'click'])
    }
    invisible(data)
}

plotCeps <- function(x, tMin, tMax, col='blue') {
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
                    yMin=0, yMax=settings$wl/2, tMin=tMin, tMax=tMax, col=col)
    }
}

plotWhistle <- function(x, tMin, tMax, plotSr, col='blue') {
    UIDs <- getWhistleData(x)$UID
    if(is.null(UIDs) ||
       length(UIDs) == 0) {
        return(NULL)
    }
    binData <- getBinaryData(x, UIDs, type='whistle')
    settings <- getPamFft(binData)
    for(i in seq_along(binData)) {
        binData[[i]]$contour <- binData[[i]]$contour * settings$sr / plotSr
        plotOneCeps(binData[[i]], hop=settings$hop, sr=settings$sr,
                    yMin=0, yMax=settings$wl/2, tMin=tMin, tMax=tMax, col=col)
    }
}

plotOneCeps <- function(x, hop, sr, yMin=0, yMax, tMin, tMax, col='blue') {
    xVals <- x$date + seq(from=0, by=hop/sr, length.out=x$nSlices)
    yVals <- x$contour

    xPlot <- scaleToOne(xVals, tMin, tMax)
    yPlot <- scaleToOne(yVals, yMin, yMax)
    lines(x=xPlot, y=yPlot, col=col, lwd=2)
}

plotClick <- function(x, tMin, tMax, yMin, yMax, col='blue') {
    clickData <- getClickData(x)[, c('UTC', 'peak', 'peakTime')]
    if(is.null(clickData) ||
       nrow(clickData) == 0) {
        return(NULL)
    }
    xVals <- clickData$UTC + clickData$peakTime
    yVals <- clickData$peak
    xPlot <- scaleToOne(xVals, tMin, tMax)
    yPlot <- scaleToOne(yVals, yMin, yMax)
    points(x=xPlot, y=yPlot, pch=1, col=col)
}

scaleToOne <- function(vals, min, max) {
    vals <- as.numeric(vals)
    min <- as.numeric(min)
    max <- as.numeric(max)
    (vals - min) / (max-min)
}
