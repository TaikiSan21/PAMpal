#' @title Calculate Average Spectra of Clicks
#'
#' @description Calculates the average spectra of all the clicks present in an
#'   event
#'
#' @param x an \linkS4class{AcousticEvent} or \linkS4class{AcousticStudy} object
#' @param evNum if \code{x} is a study, the event index number to calculate the average
#'   spectra for. Note that this is the index in the order that they appear in the
#'   \linkS4class{AcousticStudy} object, not the actual event number. Alternatively
#'   full event names can be used
#' @param calibration a calibration function to apply, if desired
#' @param wl the size of the click clips to use for calculating the spectrum. If
#'   greater than the clip present in the binary, clip will be zero padded
#' @param channel channel(s) to include in calculations. Currently does not
#'   correspond to actual channel in instrument, just the order present in the
#'   binary file
#' @param filterfrom_khz frequency in khz of highpass filter to apply, or the lower
#'   bound of a bandpass filter if \code{filterto_khz} is not \code{NULL}
#' @param filterto_khz if a bandpass filter is desired, set this as the upper bound.
#'   If only a highpass filter is desired, leave as the default \code{NULL} value.
#'   Currently only highpass and bandpass filters are supported, so if
#'   \code{filterfrom_khz} is left as zero then this parameter will have no effect
#' @param sr a sample rate to use if the sample rate present in the database needs
#'   to be overridden (typically not needed)
#' @param snr minimum signal-to-noise ratio to be included in the average, in dB. SNR is
#'   calculated as difference between the signal and noise spectra at the peak frequency
#'   of the signal. This can be inaccurate if noise is inaccurate (see \code{noise} for
#'   issues with noise calculations)
#' @param norm logical flag to normalize dB magnitude to maximum of 0
#' @param plot logical flag whether or not to plot the result. This will create two
#'   plots, the first is a concatenated spectrogram where the y-axis is
#'   frequency and the x-axis is click number. The second plot is the average
#'   spectrogram of all clicks, the y-axis is dB, x-axis
#'   is frequency.
#' @param noise logical flag to plot an average noise spectrum. This is estimated
#'   by taking a window of length \code{wl} immediately before click. Since there
#'   are only a limited number of samples saved in the Pamguard binary files, this
#'   can be very inaccurate when \code{wl} is a large proportion of the total samples
#'   saved. In these cases the noise floor will appear nearly identical to the signal,
#'   reducing \code{wl} can help get a more accurate noise floor.
#' @param \dots optional args
#'
#' @return invisibly returns a list with six items: \code{freq} - the frequency,
#'   \code{UID} - the UID of each click, \code{avgSpec} - the average spectra of the event,
#'   \code{allSpec} - the individual spectrum of each click in the event as a matrix with
#'   each spectrum in a separate column, \code{avgNoise} - the average noise spectra,
#'   \code{allNoise} - the individual noise spectrum for each click
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
#' str(avSpec$avgSpec)
#' range(avSpec$freq)
#' str(avSpec$allSpec)
#'
#' @importFrom signal hanning
#' @importFrom graphics par image axis lines legend
#' @importFrom stats fft
#'
#' @export
#'
calculateAverageSpectra <- function(x, evNum=1, calibration=NULL, wl=512,
                                    channel = 1:2, filterfrom_khz=0, filterto_khz=NULL,
                                    sr=NULL, snr=0, norm=TRUE, plot=TRUE, noise=FALSE, ...) {
    if(is.AcousticEvent(x)) {
        ev <- x
    } else if(is.AcousticStudy(x)) {
        ev <- x[evNum]
        if(is.null(sr)) {
            sr <- getSr(x)
        }
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
    # calc spectra
    specData <- lapply(binData, function(x) {
        if(is.null(x$wave)) {
            return(rep(NA, wl%/%2))
        }
        result <- 0
        hasChan <- 1:ncol(x$wave)
        useChan <- hasChan[hasChan %in% channel]
        for(c in useChan) {
            tmp <- 10^((myGram(x, channel = c, wl=wl,
                                      from=filterfrom_khz,
                                      to=filterto_khz,
                                      sr = sr, noise=FALSE, ...)$dB + calFun(freq))/20)
            if(any(is.na(tmp))) next
            result <- result + tmp
        }
        if(all(result == 0)) {
            return(rep(NA, wl%/%2))
        }
        20*log10(result / length(useChan))
    })
    specData <- specData[!sapply(specData, is.null)]
    specMat <- matrix(NA,nrow=length(specData[[1]]), ncol=length(specData))

    for(i in seq_along(specData)) {
        specMat[, i] <- specData[[i]]
    }
    snrMat <- rep(0, ncol(specMat))
    # calc noise spectra
    noiseData <- lapply(binData, function(x) {
        if(is.null(x$wave)) {
            return(rep(NA, wl%/%2))
        }
        result <- 0
        hasChan <- 1:ncol(x$wave)
        useChan <- hasChan[hasChan %in% channel]
        for(c in useChan) {
            tmp <- 10^((myGram(x, channel = c, wl=wl,
                               from=filterfrom_khz,
                               to=filterto_khz,
                               sr = sr, noise=TRUE, ...)$dB + calFun(freq))/20)
            if(any(is.na(tmp))) next
            result <- result + tmp
        }
        if(all(result == 0)) {
            return(rep(NA, wl%/%2))
        }
        20*log10(result / length(useChan))
    })
    noiseData <- noiseData[!sapply(noiseData, is.null)]
    noiseMat <- matrix(NA,nrow=length(noiseData[[1]]), ncol=length(noiseData))

    for(i in seq_along(noiseData)) {
        noiseMat[, i] <- noiseData[[i]]
    }
    # snrMat <- specMat - noiseMat
    # snrMat <- apply(snrMat, 2, function(y) {
    #     if(all(is.na(y))) return(Inf)
    #     max(y, na.rm=TRUE)
    # })

    snrVals <- vector('numeric', length = ncol(specMat))
    for(i in seq_along(snrVals)) {
        wherePeak <- which.max(specMat[, i])
        if(length(wherePeak) == 0) {
            snrVals[i] <- NA
            next
        }
        if(is.na(noiseMat[wherePeak, i])) {
            snrVals[i] <- Inf
            next
        }
        snrVals[i] <- specMat[wherePeak, i] - noiseMat[wherePeak, i]
    }
    averageNoise <- 20*log10(apply(noiseMat, 1, function(x) {
        mean(10^(x/20), na.rm=TRUE)
    }))

    averageSpec <- 20*log10(apply(specMat[, snrVals >= snr, drop=FALSE], 1, function(x) {
        mean(10^(x/20), na.rm=TRUE)
    }))

    if(norm) {
        averageNoise <- averageNoise - max(averageSpec, na.rm=TRUE)
        averageSpec <- averageSpec - max(averageSpec, na.rm=TRUE)
    }
    if(plot) {
        plotMat <- specMat[, snrVals >= snr, drop=FALSE]
        lim <- mean(plotMat, na.rm=TRUE) + c(-1,1) * 3 * sd(plotMat, na.rm=TRUE)
        plotMat[plotMat < lim[1]] <- lim[1]
        plotMat[plotMat > lim[2]] <- lim[2]
        image(t(plotMat), xaxt='n', yaxt='n', ylab='Frequency (kHz)', xlab='Click Number')
        title('Concatenated Click Spectrogram')
        xPretty <- pretty(1:ncol(plotMat), n=5)
        axis(1, at = xPretty/ncol(plotMat), labels = xPretty)
        freqPretty <- pretty(0:max(freq/1e3), n=5)
        axis(2, at = freqPretty/max(freq/1e3), labels = freqPretty)
        ylab <- ifelse(norm, 'Normalized Magnitude (dB)', 'Magnitude (dB)')
        plot(x=freq, averageSpec, type='l',
             xaxt='n', yaxt='n', ylab=ylab, xlab='Frequency (kHz)')
        if(noise) {
            lines(x=freq, averageNoise, type='l', lty=3, lwd=2)
            # legend('topright', inset=c(0, 0), xpd=TRUE, legend=c('Signal', 'Noise'), lty=1:2, bty='n')
        }
        title('Average Spectrum')
        axis(1, at = freqPretty*1e3, labels=freqPretty)
        yPretty <- pretty(range(averageSpec), n=5)
        axis(2, at=yPretty, labels=yPretty)

    }
    invisible(list(freq=freq, UID = names(specData), avgSpec=averageSpec, allSpec=specMat,
                   avgNoise=averageNoise, allNoise=noiseMat))
    # invisible(list(freq=freq, average=averageSpec, all=specMat))
}

myGram <- function(x, channel=1, wl = 512, window = TRUE, sr=NULL,
                   from=0, to=NULL, noise=FALSE) {
    wave <- x$wave[, channel]
    if(is.null(sr)) {
        sr <- x$sr
    }
    if(from > 0) {
        # kinda janky because NULL * 1e3 is not NULL anymore, its numeric(0)
        if(!is.null(to)) {
            to_hz <- to * 1e3
        } else {
            to_hz <- NULL
        }
        wave <- bwfilter(wave, f=sr, n=4, from=from*1e3, to = to_hz, output='sample')
    }
    wave <- clipAroundPeak(wave, wl, noise=noise)

    FUN <- function(x) {
        result <- Mod(fft(x))
        20*log10(result)
    }

    y <- (1:(wl)) / wl * sr

    if(window) {
        wave <- wave * hanning(length(wave)) / mean(hanning(length(wave)))
    }
    dB <- FUN(wave)[1:(wl)]
    isInf <- is.infinite(dB)
    if(all(isInf)) {
        return(list(dB = rep(NA, wl%/%2), freq=y[1:(wl%/%2)]))
    }
    if(any(isInf)) {
        dB[dB == Inf] <- max(dB[!isInf], na.rm=TRUE)
        dB[dB == -Inf] <- min(dB[!isInf], na.rm=TRUE)
    }
    list(dB = dB[1:(wl%/%2)], freq=y[1:(wl%/%2)])
}

getSr <- function(x) {
    clickFuns <- pps(x)@functions$ClickDetector
    srSettings <- unique(unlist(lapply(clickFuns, function(c) {
        formals(c)[['sr_hz']]
    })))
    if(length(srSettings) == 1 &&
       is.numeric(srSettings)) {
        return(srSettings)
    }
    NULL
}
