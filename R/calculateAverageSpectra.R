#' @title Calculate Average Spectra of Clicks
#'
#' @description Calculates the average spectra of all the clicks present in an
#'   event
#'
#' @param x an \linkS4class{AcousticEvent} or \linkS4class{AcousticStudy} object
#' @param evNum if \code{x} is a study, the event number to calculate the average
#'   spectra for
#' @param calibration a calibration function to apply, if desired
#' @param wl the size of the click clips to use for calculating the spectra. If
#'   greater than the clip present in the binary, clip will be zero padded
#' @param sr a sample rate to use if the sample rate present in the database needs
#'   to be overridden (typically only needed if a decimator was used)
#' @param plot logical flag whether or not to plot the result. The plot will be a
#'   two panel plot, the top is a concatenated spectrogram where the y-axis is
#'   frequency and the x-axis is click number. The bottom plot is the average
#'   spectrogram of all clicks, the y-axis is normalized magnitude (dB values
#'   for each click are normalized between 0 and 1 before averaging), x-axis
#'   is frequency.
#'
#' @return invisibly returns a list with three items: \code{freq}, the frequency,
#'   \code{average}, the average spectra of the event, and \code{all}, the individual
#'   spectra of each click in the event as a matrix.
#'
#' @author Taiki Sakai \email{taiki.sakai@@noaa.gov}
#'
#' @importFrom signal hanning
#' @importFrom graphics par image axis
#' @importFrom stats fft
#' @export
#'
calculateAverageSpectra <- function(x, evNum=1, calibration=NULL, wl=4096, sr=NULL, plot=TRUE) {
    if(is.AcousticEvent(x)) {
        ev <- x
    } else if(is.AcousticStudy(x)) {
        ev <- events(x)[[evNum]]
    }
    clicks <- getDetectorData(ev)$click
    binData <- getBinaryData(ev, clicks$UID, quiet = TRUE)
    if(length(binData) == 0) {
        stop('Not able to find any data for this event.')
    }
    if(is.null(calibration)) {
        calFun <- function(x) 0
    } else if(is.function(calibration)) {
        calFun <- calibration
    } else if(is.character(calibration)) {
        calFun <- findCalibration(calibration)
    }
    freq <- myGram(binData[[1]], channel=1, wl=wl, sr=sr)$freq
    specData <- lapply(binData, function(x) {
        myGram(x, channel = 1, wl=wl)$dB + calFun(freq)
    })
    specMat <- matrix(NA,nrow=length(specData[[1]]), ncol=length(specData))
    for(i in seq_along(specData)) {
        specMat[, i] <- (specData[[i]] - min(specData[[i]]))/diff(range(specData[[i]]))
    }
    averageSpec <- apply(specMat, 1, mean)
    if(plot) {
        par(mfrow=c(2,1))
        image(t(specMat), xaxt='n', yaxt='n', ylab='Frequency (kHz)', xlab='Click Number')
        xPretty <- pretty(1:length(specData), n=5)
        axis(1, at = xPretty/length(specData), labels = xPretty)
        yPretty <- pretty(0:max(freq/1e3), n=5)
        axis(2, at = yPretty/max(freq/1e3), labels = yPretty)
        plot(x=freq, averageSpec, type='l',
             xaxt='n', yaxt='n', ylab='Normalized Magnitude', xlab='Frequency (kHz)')
        axis(1, at = yPretty*1e3, labels=yPretty)
        yPretty <- pretty(0:1, n=5)
        axis(2, at=yPretty, labels=yPretty)
    }
    invisible(list(freq=freq, average=averageSpec, all=specMat))
}

myGram <- function(x, channel=1, wl = 2048, window = T, sr=NULL) {
    wave <- x$wave[, channel]
    if(is.null(sr)) {
        sr <- x$sr
    }
    wave <- clipAroundPeak(wave, wl)
    FUN <- function(x) {
        20*log10(abs(fft(x)))
    }

    y <- (1:(wl)) / wl * sr

    if(window) wave <- wave * hanning(length(wave)) / mean(hanning(length(wave)))
    dB <- FUN(wave)[1:(wl)]
    list(dB = dB[1:(wl%/%2)], freq=y[1:(wl%/%2)])
}
