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
#' @param norm logical flag to normalize magnitudes to 0-1 range
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
#' @examples
#'
#' data(exStudy)
#' # need to update binary file locations to users PAMpal installation
#' binUpd <- system.file('extdata', 'Binaries', package='PAMpal')
#' dbUpd <- system.file('extdata', package='PAMpal')
#' exStudy <- updateFiles(exStudy, bin = binUpd, db=dbUpd)
#' avSpec <- calculateAverageSpectra(exStudy)
#' str(avSpec$average)
#' range(avSpec$freq)
#' str(avSpec$all)
#'
#' @importFrom signal hanning
#' @importFrom graphics par image axis
#' @importFrom stats fft
#' @export
#'
calculateAverageSpectra <- function(x, evNum=1, calibration=NULL, wl=1024,
                                    sr=NULL, norm=TRUE, plot=TRUE) {
    if(is.AcousticEvent(x)) {
        ev <- x
    } else if(is.AcousticStudy(x)) {
        ev <- events(x)[[evNum]]
        if(is.null(calibration) &&
           length(x@pps@calibration$ClickDetector) == 1) {
            calibration <- x@pps@calibration$ClickDetector[[1]]
        }
    }
    clicks <- getDetectorData(ev)$click
    binData <- getBinaryData(ev, clicks$UID, quiet = TRUE)
    if(length(binData) == 0) {
        stop('Not able to find any data for this event.')
    }
    if(is.null(calibration) ||
       length(calibration) == 0) {
        calFun <- function(x) 0
    } else if(is.function(calibration)) {
        calFun <- calibration
    } else if(is.character(calibration)) {
        calFun <- findCalibration(calibration)
    } else {
        calFun <- function(x) 0
    }
    freq <- myGram(binData[[1]], channel=1, wl=wl, sr=sr)$freq
    specData <- lapply(binData, function(x) {
        if(is.null(x$wave)) {
            return(NULL)
        }
        result <- 0
        for(c in 1:ncol(x$wave)) {
            result <- result + myGram(x, channel = c, wl=wl)$dB + calFun(freq)
        }
        result / ncol(x$wave)
    })
    specData <- specData[!sapply(specData, is.null)]
    specMat <- matrix(NA,nrow=length(specData[[1]]), ncol=length(specData))

    for(i in seq_along(specData)) {
        if(norm) {
            specMat[, i] <- (specData[[i]] - min(specData[[i]]))/diff(range(specData[[i]]))
        } else {
            specMat[, i] <- specData[[i]]
        }
    }

    averageSpec <- apply(specMat, 1, mean)
    if(plot) {
        oldMf <- par()$mfrow
        par(mfrow=c(2,1))
        image(t(specMat), xaxt='n', yaxt='n', ylab='Frequency (kHz)', xlab='Click Number')
        xPretty <- pretty(1:length(specData), n=5)
        axis(1, at = xPretty/length(specData), labels = xPretty)
        freqPretty <- pretty(0:max(freq/1e3), n=5)
        axis(2, at = freqPretty/max(freq/1e3), labels = freqPretty)
        ylab <- ifelse(norm, 'Normalized Magnitude', 'Magnitude (dB)')
        plot(x=freq, averageSpec, type='l',
             xaxt='n', yaxt='n', ylab=ylab, xlab='Frequency (kHz)')
        axis(1, at = freqPretty*1e3, labels=freqPretty)
        yPretty <- pretty(range(averageSpec), n=5)
        axis(2, at=yPretty, labels=yPretty)
        par(mfrow = oldMf)
    }
    invisible(list(freq=freq, average=averageSpec, all=specMat))
}

myGram <- function(x, channel=1, wl = 512, window = TRUE, sr=NULL) {
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
    isInf <- is.infinite(dB)
    if(any(isInf)) {
        dB[dB == Inf] <- max(dB[!isInf])
        dB[dB == -Inf] <- min(dB[!isInf])
    }
    list(dB = dB[1:(wl%/%2)], freq=y[1:(wl%/%2)])
}
