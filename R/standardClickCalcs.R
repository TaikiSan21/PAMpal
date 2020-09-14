#' @title Calculate a Set of Measurements for Clicks
#'
#' @description Calculate a set of "standard" measurements for odontocete clicks.
#'   Most of calculations following approach of Baumann-Pickering / Soldevilla
#'
#' @param data a list that must have 'wave' containing the wave form as a
#'   matrix with a separate column for each channel, and 'sr' the
#'   sample rate of the data. Data can also be a \code{Wave} class
#'   object, like one created by \code{\link[tuneR]{Wave}}.
#' @param sr_hz either \code{'auto'} (default) or the numeric value of the sample
#'   rate in hertz. If \code{'auto'}, the sample rate will be read from the
#'   'sr' of \code{data}
#' @param calibration a calibration function to apply to the spectrum, must be
#'   a gam. If NULL no calibration will be applied (not recommended).
#' @param filterfrom_khz frequency in khz of highpass filter to apply, or the lower
#'   bound of a bandpass filter if \code{filterto_khz} is not \code{NULL}
#' @param filterto_khz if a bandpass filter is desired, set this as the upper bound.
#'   If only a highpass filter is desired, leave as the default \code{NULL} value
#' @param winLen_sec length in seconds of fft window. The click wave is first
#'   shortened to this number of samples around the peak of the wave,
#'   removing a lot of the noise around the click. Following approach of
#'   JB/EG/MS.
#'
#' @return A data frame with one row for each channel of click waveform.
#'   Calculates approximate noise level and click duration from the
#'   TK energy (Soldevilla JASA17), up to 3 highest peak frequencies and
#'   the 'troughs' between them (see \code{\link[PAMmisc]{peakTrough}}), and the 3
#'   and 10dB bandwidth levels and 'Q' value (see \code{\link[seewave]{Q}}).
#'
#' @author Taiki Sakai \email{taiki.sakai@@noaa.gov}
#'
#' @importFrom seewave bwfilter TKEO spec
#' @importFrom dplyr bind_rows
#' @importFrom PAMmisc peakTrough
#' @importFrom stats median quantile
#' @importFrom tuneR WaveMC
#' @export
#'
standardClickCalcs <- function(data, sr_hz='auto', calibration=NULL, filterfrom_khz=10, filterto_khz=NULL, winLen_sec=.0025) {
    # SLOWEST PART BY FAR is bwfilter
    result <- list()
    paramNames <- c('Channel', 'noiseLevel', 'duration', 'peakTime', 'peak', 'peak2', 'peak3', 'trough',
                    'trough2', 'peakToPeak2', 'peakToPeak3', 'peak2ToPeak3', 'dBPP', 'Q_10dB',
                    'PeakHz_10dB', 'fmin_10dB', 'fmax_10dB', 'BW_10dB', 'centerHz_10dB',
                    'Q_3dB', 'PeakHz_3dB', 'fmin_3dB', 'fmax_3dB', 'BW_3dB', 'centerHz_3dB')
    # Do for each channel
    if(inherits(data, 'Wave')) {
        data <- WaveMC(data)
        data <- list(wave = data@.Data, sr = data@samp.rate)
    }
    if(!is.matrix(data$wave)) {
        data$wave <- matrix(data$wave, ncol=1)
    }
    for(chan in 1:ncol(data$wave)) {
        # We store results in 'thisDf', note channels start at 1 not 0
        thisDf <- list(Channel = chan)
        thisWave <- data$wave[,chan]
        if(all(thisWave == 0)) {
            # cant return NULL in this case - if other functions compute something useful
            # we need to be able to bind_cols which requires same number of rows
            blanks  <- data.frame(matrix(NA, nrow=1, ncol=length(paramNames)))
            colnames(blanks) <- paramNames
            blanks[['Channel']] <- chan
            result[[chan]] <- blanks
            next
        }
        if(sr_hz == 'auto') {
            sr <- data$sr
        } else {
            sr <- sr_hz
        }

        if(filterfrom_khz > 0) {
            # kinda janky because NULL * 1e3 is not NULL anymore, its numeric(0)
            if(!is.null(filterto_khz)) {
                to_hz <- filterto_khz * 1e3
            } else {
                to_hz <- NULL
            }
            thisWave <- bwfilter(thisWave, f=sr, n=4, from=filterfrom_khz*1e3, to = to_hz, output='sample')
        }
        # -1 here because time after start time
        peakTime <- (which.max(abs(thisWave)) - 1) / sr

        # 2.5ms window size - following Soldevilla paper JASA17
        fftSize <- round(sr * winLen_sec, 0)
        fftSize <- fftSize + (fftSize %% 2)

        # Shortening click wave by getting numpoints=fftsize around peak of wav. Length fft+1
        # This removes a lot of the noise samples Pamguard takes - following what Jay/EG did
        # wavPeak <- which.max(thisWave)

        # if(length(thisWave) <= fftSize) {
        #     # do nothing
        # } else if(wavPeak > fftSize/2) {
        #     if((wavPeak + fftSize/2) <= length(thisWave)) {
        #         thisWave <- thisWave[(wavPeak - fftSize/2):(wavPeak + fftSize/2)]
        #     } else {
        #         thisWave <- thisWave[(length(thisWave)-fftSize):length(thisWave)]
        #     }
        # } else {
        #     thisWave <- thisWave[1:(fftSize+1)]
        # }
        thisWave <- clipAroundPeak(thisWave, fftSize)

        # Do any calculations you want. Here just getting peak frequency.
        # TKEO - skip 1st .001s to avoid startup artifacts, energy is 2nd col
        # This .001s bit doesnt work for really short samples...
        # UPDATE 8-24-2020 THERES NO WAY THE STARTUP ARTIFACT IS RELEVANT
        thisTk <- TKEO(thisWave, f=sr, M=1,plot=F)
        tkEnergy <- thisTk[1:length(thisWave),2]
        tkDb <- 10*log10(tkEnergy-min(tkEnergy, na.rm=TRUE))
        tkDb <- tkDb - max(tkDb, na.rm=TRUE)
        tkDb[!is.finite(tkDb)] <- NA

        noiseLevel <- median(tkDb, na.rm=TRUE)
        if(is.na(noiseLevel)) {
            noiseLevel <- 0
        }
        thisDf$noiseLevel <- noiseLevel

        # duration defined as 100 time above 40% TKE threshold
        noiseThresh <- quantile(thisTk[,2], probs=.4, na.rm=TRUE)*100
        dur <- subset(thisTk, thisTk[,2] >= noiseThresh)
        if(length(dur)==0) {
            dur <- 0
        } else {
            dur <- 1e6*(max(dur[,1])-min(dur[,1]))
        }
        thisDf$duration <- dur
        thisDf$peakTime <- peakTime
        thisSpec <- spec(thisWave, f=sr, wl=fftSize, norm=FALSE, correction='amplitude', plot=FALSE)
        if(any(is.nan(thisSpec[,2]))) {
            blanks  <- data.frame(matrix(NA, nrow=1, ncol=length(paramNames)))
            colnames(blanks) <- paramNames
            blanks[['Channel']] <- chan
            result[[chan]] <- blanks
            next
        }
        relDb <- 20*log10(thisSpec[,2])
        # only put this in an IF so its easy to have a breakpoint
        if(any(!is.finite(relDb))) {
            relDb[!is.finite(relDb)] <- NA
        }
        # This is for integer overflow spec jankiness fixing
        freq <- seq(from=0, by = thisSpec[2,1] - thisSpec[1,1], length.out = nrow(thisSpec))
        # Calibration - I don't have a standardized way of doing this yet
        # newClick <- data.frame(Freq=thisSpec[,1]*1e3, Sens = relDb)


        if(!is.null(calibration)) {
            #DO CA
            if(is.function(calibration)) {
                calFun <- calibration
            } else if(is.character(calibration)) {
                calFun <- findCalibration(calibration)
            }
            relDb <- relDb + calFun(freq * 1e3)
        }
        # } else {
        #     # if no cal, just use original relDb
        #     clickSens <- relDb - max(relDb, na.rm=TRUE)
        # }
        calibratedClick <- cbind(freq, relDb)

        peakData <- peakTrough(calibratedClick)
        thisDf <- c(thisDf, peakData)

        # peak-to-peak calcuation from anne s
        dBPP <- 20*log10(max(thisWave) - min(thisWave))
        if(!is.null(calibration)) {
            dBPP <- dBPP + calFun(peakData$peak * 1e3)
        }

        thisDf$dBPP <- dBPP
        # Finding 10/3 dB bandwidth - modified 'Q' function from seewave package
        dbBW10 <- Qfast(calibratedClick, f=sr, level=-10, plot=FALSE)
        names(dbBW10) <- c('Q_10dB', 'PeakHz_10dB', 'fmin_10dB', 'fmax_10dB', 'BW_10dB')
        dbBW10$centerHz_10dB <- dbBW10$fmax_10dB - (dbBW10$BW_10dB/2)

        dbBW3 <- Qfast(calibratedClick, f=sr, level=-3, plot=FALSE)
        names(dbBW3) <- c('Q_3dB', 'PeakHz_3dB', 'fmin_3dB', 'fmax_3dB', 'BW_3dB')
        dbBW3$centerHz_3dB <- dbBW3$fmax_3dB - (dbBW3$BW_3dB/2)

        thisDf <- c(thisDf, dbBW10, dbBW3)

        result[[chan]] <- thisDf
    }
    # Combine calcs for all channels
    result <- bind_rows(result)
    result$Channel <- as.character(result$Channel)
    result
}

#' @importFrom stats approx
#'
Qfast <- function(spec,
                  f = NULL,
                  level = -10,
                  mel = FALSE,
                  plot = FALSE,
                  colval = "red",
                  cexval = 1,
                  fontval = 1,
                  flab = NULL,
                  alab = "Relative amplitude (dB)",
                  type = "l", ...) {

    if (!is.null(f) & mel) {
        f <- 2 * mel(f/2)
    }
    if (is.null(f)) {
        if (is.vector(spec))
            stop("'f' is missing")
        else if (is.matrix(spec))
            f <- spec[nrow(spec), 1] * 2000 * nrow(spec)/(nrow(spec) - 1)
    }
    if (is.matrix(spec)) {
        range <- c(spec[1,1], spec[nrow(spec),1]) ### THERES NO RANGE IF SPEC IS VEC
        # spectest <- spec
        spec <- spec[, 2]
    }
    specMax <- which.max(spec)
    if(length(specMax)==0) {
        return(list(Q = 0, dfreq = 0, fmin = 0, fmax = 0, bdw = 0))
    }
    # if (spec[specMax] == 1)
    #     stop("data must be in dB")
    # if (specMax == 1)
    #     stop("maximal peak cannot be the first value of the spectrum")

    n1 <- length(spec)
    level2 <- spec[specMax] + level
    f0 <- specMax
    f0khz <- (((f0-1)/(n1-1)))*(range[2]-range[1]) + range[1] # These and others below need -1s to properly scale
    specA <- spec[1:f0]
    specB <- spec[f0:length(spec)]
    negA <- which(specA <= level2)
    if(length(negA) == 0) {
        fA <- 1
    } else {
        firstNegA <- max(which(specA <= level2)) # btwn this and next
        # fA <- approx(x=spec[firstNegA:(firstNegA+1)], y=firstNegA:(firstNegA+1), xout=level2)$y
        fA <- oneInterp(x=spec[firstNegA:(firstNegA+1)], y=firstNegA:(firstNegA+1), xout=level2)
    }
    fAkhz <- ((fA-1)/(n1-1)) * (range[2]-range[1]) + range[1]
    nA <- length(specA)
    negB <- which(specB <= level2)
    if(length(negB) == 0) {
        fB <- length(spec)
    } else {
        firstNegB <- min(negB) + (nA - 1) # btwn this and prev
        # fB <- approx(x=spec[(firstNegB-1):firstNegB], y=(firstNegB-1):firstNegB, xout=level2)$y
        fB <- oneInterp(x=spec[(firstNegB-1):firstNegB], y=(firstNegB-1):firstNegB, xout=level2)
    }
    fBkhz <- ((fB-1)/(n1-1)) * (range[2]-range[1]) + range[1]
    Q <- f0khz/(fBkhz-fAkhz)
    results <- list(Q = Q, dfreq = f0khz, fmin = fAkhz, fmax = fBkhz,
                    bdw = fBkhz - fAkhz)

    # Temp fix on missing, if any are borked then Q is borked so go check
    if(length(Q) == 0) {
        results <- lapply(results, function(x) ifelse(length(x)==0, 0, x))
    }
    # if (plot) {
    #     if (is.null(flab)) {
    #         if (mel)
    #             flab <- "Frequency (kmel)"
    #         else flab <- "Frequency (kHz)"
    #     }
    #     x <- seq(range[1], range[2], length.out = n1)
    #     plot(x = x, y = spec, xlab = flab, ylab = alab, type = type,
    #          ...)
    #     arrows(fAkhz, level2, fBkhz, level2, length = 0.1, col = colval,
    #            code = 3, angle = 15)
    #     text(paste("Q =", as.character(round(Q, 2))), x = fBkhz,
    #          y = level2, pos = 4, col = colval, cex = cexval,
    #          font = fontval)
    #     invisible(results)
    # }
    return(results)
}

oneInterp <- function(x, y, xout) {
    (y[1]*(x[2]-xout) - y[2]*(x[1]-xout)) / (x[2]-x[1])
}
