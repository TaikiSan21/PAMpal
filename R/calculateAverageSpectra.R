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
#'   spectrogram of all clicks, the y-axis is dB, x-axis is frequency. Can be a
#'   vector of length two to create only one of the two plots
#' @param noise logical flag to plot an average noise spectrum. This is estimated
#'   by taking a window of length \code{wl} immediately before click. Since there
#'   are only a limited number of samples saved in the Pamguard binary files, this
#'   can be very inaccurate when \code{wl} is a large proportion of the total samples
#'   saved. In these cases the noise floor will appear nearly identical to the signal,
#'   reducing \code{wl} can help get a more accurate noise floor.
#' @param decimate integer factor to reduce sample rate by
#' @param sort logical flag to sort concatenated spectrogram by peak frequency
#' @param mode one of \code{'spec'} or \code{'ceps'} to plot the spectrum or cepstrum
#' @param title replacement titles for plots, can be length vector of length two to provide
#'   separate titles
#' @param ylim optional y limits for mean spectra plot
#' @param flim optional frequency limits for both plots
#' @param cmap colors to use for concatenated click spectrogram, either a palette
#'   function or vector of colors
#' @param brightness value from -255 to 255, positive values increase brightness,
#'   negative values decrease brightness of concatenated spectrogram image
#' @param contrast value from -255 to 255, positive values increase contrast, negative
#'   values decrease contrast of concatenated spectrogram image
#' @param q lower and upper quantiles to remove for scaling concatenated spectrogram.
#'   Or if a single value, then quantiles \code{q} and \code{1-q} will be used. Ex.
#'   if \code{q=.01}, then the bottom 1% and top 1% of values are truncated before
#'   plotting the image. This is done purely for cosmetic reasons, no output data is
#'   affected
#' @param showBreaks logical flag to show lines separating events when plotting
#'   multiple events
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
#' @importFrom signal hanning decimate
#' @importFrom graphics par image axis lines legend
#' @importFrom stats fft
#' @importFrom grDevices hcl.colors
#'
#' @export
#'
calculateAverageSpectra <- function(x, evNum=1, calibration=NULL, wl=512,
                                    channel = 1:2, filterfrom_khz=0, filterto_khz=NULL,
                                    sr=NULL, snr=0, norm=TRUE, plot=TRUE, noise=FALSE,
                                    decimate=1,
                                    sort=FALSE, mode='spec',
                                    title=NULL,
                                    ylim=NULL,
                                    flim=NULL,
                                    cmap=hcl.colors(30, "YlOrRd", rev = TRUE),
                                    brightness=0,
                                    contrast=0,
                                    q=.01,
                                    showBreaks=TRUE,
                                    ...) {
    if(is.AcousticEvent(x)) {
        ev <- x
    } else if(is.AcousticStudy(x)) {
        if(is.numeric(evNum)) {
            evNum <- evNum[evNum <= length(events(x))]
        }
        ev <- x[evNum]
        if(is.null(sr)) {
            sr <- getSr(x, type='click')
        }
        if(is.null(calibration) &&
           length(x@pps@calibration$ClickDetector) == 1) {
            calibration <- x@pps@calibration$ClickDetector[[1]]
        }
    }
    clickData <- getClickData(ev)[c('eventId', 'UID', 'UTC')]
    if(is.null(clickData)) {
        stop('No clicks in this event')
        return(NULL)
    }
    clickData <- arrange(clickData, .data$UTC)[c('eventId', 'UID')]
    clickData <- distinct(clickData)
    clickUID <- clickData$UID
    # reordering because they get put in UID order
    binData <- lapply(unname(split(clickData, clickData$eventId)), function(c) {
        bin <- getBinaryData(ev[[c$eventId[1]]], c$UID, type='click', quiet=TRUE)
        bin[c$UID]
    })
    binData <- unlist(binData, recursive = FALSE)

    # binData <- getBinaryData(ev, clickUID, type='click', quiet = TRUE)[clickUID]
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
    freq <- myGram(binData[[1]], channel=1, wl=wl, sr=sr, mode=mode,
                   decimate=decimate)[, 1]

    specData <- binToSpecMat(binData, channel=channel, wl=wl, sr=sr, freq=freq,
                             calFun=calFun, mode=mode, noise=FALSE, decimate=decimate,
                             filterfrom_khz=filterfrom_khz, filterto_khz = filterto_khz, ...)
    specData <- specData[!sapply(specData, is.null)]
    specMat <- matrix(NA,nrow=length(specData[[1]]), ncol=length(specData))

    for(i in seq_along(specData)) {
        specMat[, i] <- specData[[i]]
    }
    snrMat <- rep(0, ncol(specMat))

    noiseData <- binToSpecMat(binData, channel=channel, wl=wl, sr=sr, freq=freq,
                              calFun=calFun, mode=mode, noise=TRUE, decimate=decimate,
                              filterfrom_khz=filterfrom_khz, filterto_khz = filterto_khz, ...)
    noiseData <- noiseData[!sapply(noiseData, is.null)]
    noiseMat <- matrix(NA,nrow=length(noiseData[[1]]), ncol=length(noiseData))

    for(i in seq_along(noiseData)) {
        noiseMat[, i] <- noiseData[[i]]
    }
    # averageNoise <- 20*log10(apply(noiseMat, 1, function(x) {
    #     mean(10^(x/20), na.rm=TRUE)
    # }))
    averageNoise <- doLogAvg(noiseMat, log=mode=='spec')
    if(snr > 0) {
        snrVals <- vector('numeric', length = ncol(specMat))
        for(i in seq_along(snrVals)) {
            minX <- switch(mode,
                           'spec' = 1,
                           # 'ceps' = sum(freq < .0004)
                           'ceps' = min(which(specMat[-(1:2), i] < 3*median(specMat[-(1:2), i]))) + 2
            )
            wherePeak <- which.max(specMat[minX:nrow(specMat), i]) + minX - 1
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

        if(!any(snrVals >= snr)) {
            if(is.AcousticEvent(x)) {
                evId <- id(x)
            }
            if(is.AcousticStudy(x)) {
                evId <- id(x[[evNum]])
            }
            warning('No clicks above SNR threshold for event ', evId, call.=FALSE)
            return(NULL)
        }
        snrKeep <- snrVals >= snr
    } else {
        snrKeep <- rep(TRUE, ncol(specMat))
    }
    # averageSpec <- 20*log10(apply(specMat[, snrKeep, drop=FALSE], 1, function(x) {
    #     mean(10^(x/20), na.rm=TRUE)
    # }))
    averageSpec <- doLogAvg(specMat[, snrKeep, drop=FALSE], log=mode == 'spec')
    if(norm) {
        maxVal <- max(averageSpec, na.rm=TRUE)
        averageNoise <- averageNoise - maxVal
        averageSpec <- averageSpec - maxVal
    }
    if(length(plot) == 1) {
        plot <- rep(plot, 2)
    }
    if(any(plot)) {
        switch(mode,
               'spec' = {
                   unitLab <- 'Frequency (kHz)'
                   unitFactor <- 1/1e3
                   title1 <- 'Concatenated Click Spectrogram'
                   title2 <- 'Average Spectrum'
                   maxZ <- 3
               },
               'ceps' = {
                   unitLab <- 'ICI (ms)'
                   unitFactor <- 1e3
                   title1 <- 'Concatenated Click Cepstrogram'
                   title2 <- 'Average Cepstrum'
                   maxZ <- 5
               })
        if(!is.null(title)) {
            if(length(title) == 1) {
                title <- rep(title, 2)
            }
            title1 <- title[1]
            title2 <- title[2]
        }
        # change y axis
        if(is.null(flim)) {
            plotFreq <- rep(TRUE, length(freq))
        } else {
            plotFreq <- freq >= flim[1] & freq <= flim[2]
        }
        freqPretty <- pretty(seq(from=min(freq[plotFreq])*unitFactor, to=max(freq[plotFreq])*unitFactor, length.out=5), n=5)

        # Concat Spec Plot
        if(plot[1]) {
            plotMat <- specMat[plotFreq, snrKeep, drop=FALSE]
            # if(mode == 'ceps') {
            #     skips <- 1:(0.01 * nrow(plotMat))
            #     plotMat <- plotMat[-skips, ,drop=FALSE]
            # }
            # lim <- mean(plotMat, na.rm=TRUE) + c(-1,1) * maxZ * sd(plotMat, na.rm=TRUE)
            goodFreq <- freq[plotFreq] >= filterfrom_khz * 1e3
            if(!is.null(filterto_khz)) {
                goodFreq <- (freq[plotFreq] <= filterto_khz * 1e3) &
                    goodFreq
            }
            plotMat[plotMat < min(plotMat[goodFreq, ])] <- min(plotMat[goodFreq,])
            if(length(q) == 1) {
                q <- c(q, 1-q)
            }
            lim <- quantile(plotMat, q)
            plotMat[plotMat < lim[1]] <- lim[1]
            plotMat[plotMat > lim[2]] <- lim[2]
            plotMat <- plotMat / diff(range(plotMat, na.rm=TRUE)) * 255
            plotMat <- plotMat - min(plotMat, na.rm=TRUE)
            plotMat[plotMat > 255] <- 255
            plotMat[plotMat < 0] <- 0
            if(sort) {
                sortIx <- sort(apply(plotMat, 2, which.max), index.return=TRUE)$ix
                plotMat <- plotMat[, sortIx, drop=FALSE]
            }
            plotMat <- doContrast(doBrightness(plotMat, brightness), contrast)
            if(is.function(cmap)) {
                cmap <- cmap(30)
            }
            image(x=1:ncol(plotMat),
                  y=seq(from=min(freq[plotFreq]*unitFactor), to=max(freq[plotFreq]*unitFactor), length.out=nrow(plotMat)),
                  z=t(plotMat), xaxt='n', yaxt='n', ylab=unitLab, xlab='Click Number', useRaster=TRUE,
                  col = cmap,
                  breaks=seq(from=0, to=255, length.out=31))
            title(title1)
            xPretty <- pretty(1:ncol(plotMat), n=5)
            axis(1, at = xPretty, labels = xPretty)
            axis(2, at = freqPretty, labels = freqPretty)
            if(isTRUE(showBreaks) & isFALSE(sort) & length(evNum) > 1) {
                evBreaks <- table(clickData$eventId)[sapply(events(ev), id)]
                evStarts <- cumsum(evBreaks)
                for(b in 1:(length(evStarts)-1)) {
                    lines(x=rep(evStarts[b]+0.5, 2),
                          y=c(0, max(freq[plotFreq]*unitFactor)),
                          lwd=2)
                }
            }
        }
        # Mean Spec Plot
        if(plot[2]) {
            ylab <- ifelse(norm, 'Normalized Magnitude (dB)', 'Magnitude (dB)')
            plot(x=freq, averageSpec, type='l',
                 xaxt='n', yaxt='n', ylab=ylab, xlab=unitLab, lwd=1, ylim=ylim, xlim=range(freq[plotFreq]))
            if(noise) {
                lines(x=freq, averageNoise, type='l', lty=3, lwd=2)
                # legend('topright', inset=c(0, 0), xpd=TRUE, legend=c('Signal', 'Noise'), lty=1:2, bty='n')
            }

            title(title2)
            axis(1, at = freqPretty/unitFactor, labels=freqPretty)
            yRange <- range(averageSpec)
            if(!is.null(ylim)) {
                if(ylim[1] > ylim[2]) {
                    ylim <- rev(ylim)
                }
                yRange[1] <- ifelse(yRange[1] < ylim[1], ylim[1], yRange[1])
                yRange[2] <- ifelse(yRange[2] < ylim[2], ylim[2], yRange[2])
            }
            yPretty <- pretty(yRange, n=5)
            axis(2, at=yPretty, labels=yPretty)
        }

    }
    invisible(list(freq=freq, UID = names(specData), avgSpec=averageSpec, allSpec=specMat,
                   avgNoise=averageNoise, allNoise=noiseMat))
}

doLogAvg <- function(x, log=TRUE) {
    if(is.null(dim(x)) ||
       dim(x)[2] == 1) {
        # x <- matrix(x, ncol=1)
        return(x)
    }
    if(dim(x)[2] == 2) {
        if(isTRUE(log)) {
            x <- 10^(x/20)
        }
        outs <- (x[, 1] + x[, 2]) / 2
        if(isTRUE(log)) {
            outs <- 20*log10(outs)
        }
        return(outs)
    }
    if(isTRUE(log)) {
        return(20*log10(apply(x, 1, function(y) {
            mean(10^(y/20), na.rm=TRUE)
        })))
    }
    apply(x, 1, function(y) {
        mean(y, na.rm=TRUE)
    })
}
