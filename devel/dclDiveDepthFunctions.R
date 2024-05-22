library(tuneR)
library(fftw)
library(signal)
library(geosphere)
library(ggplot2)
library(PAMpal)
library(shiny)

# peak finder in R
# Copyright Nathanael C. Yoder 2015 (nyoder@gmail.com) adapted MATLAB code
# PAMmisc
peakfinder <- function(x0, thresh=NULL, extrema=1, plot=FALSE) {
    if(!is.vector(x0)) {
        error('PEAKFINDER:Input\nThe input data must be a vector')
    }
    if(is.null(x0) || length(x0) == 0) {
        return(list(peakLoc=NULL, peakMag=NULL))
    }
    if(is.complex(x0)) {
        warning('PEAKFINDER:NotReal\nAbsolute value of data will be used')
        x0 <- abs(x0)
    }
    if(is.null(thresh)) {
        thresh <- (max(x0) - min(x0)) / 4
    } else if(!is.numeric(thresh) || is.complex(thresh)) {
        thresh <- (max(x0) - min(x0)) / 4
        warning('PEAKFINDER:VectorThresh\n',
                'The threshold must be a real scalar. A threshold of ',
                round(thresh, 4), ' will be used.')
    } else if(length(thresh) > 1) {
        warning('PEAKFINDER:VectorThresh\n',
                'The threshold must be a scalar. The first threshold in the vector will be used.')
        thresh <- thresh[1]
    }
    extrema <- sign(extrema) # should only be 1 or -1
    if(extrema == 0) {
        stop('PEAKFINDER:ZeroMaxima\n',
             'Either 1 (for maxima) or -1 (for minima) must be input for extrema')
    }
    x0 <- x0 * extrema # make it so always finding max
    dx0 <- diff(x0) # find derivative
    eps <- 1e-16 
    dx0[dx0 == 0] <- -1*eps # this is so we findthe first of repeated values
    ind <- which(dx0[1:(length(dx0)-1)] * dx0[2:length(dx0)] < 0) + 1 # find where sign changes
    
    # include endpoints in potential peaks and valleys
    x <- c(x0[1], x0[ind], x0[length(x0)])
    ind <- c(1, ind, length(x0))
    
    # x only has the peaks, valleys, and endpoints
    len <- length(x)
    minMag <- min(x)
    
    if(len > 2) { # function with peaks and valleys
        # set initial parameters for loop
        tempMag <- minMag
        foundPeak <- FALSE
        leftMin <- minMag
        
        # deal with first point a little differently since tacked it on
        # calculate the sign of the derivative isnce we tacked the first point
        # on it does not necessarily alternate like the rest.
        signDx <- sign(diff(x[1:3]))
        if(signDx[1] <= 0) { # first point is larger or equal to the second
            ii <- 0
            if(signDx[1] == signDx[2]) { # want alternating signs
                x <- x[-2]
                ind <- ind[-2]
                len <- len-1
            }
        } else { # first point smaller than second
            ii <- 1
            if(signDx[1] == signDx[2]) { # want alternating signs
                x <- x[-1]
                ind <- ind[-1]
                len <- len-1
            }
        }
        
        # prealloc max number of maxima
        maxPeaks <- ceiling(len/2)
        peakLoc <- rep(0, length(maxPeaks))
        peakMag <- peakLoc
        cInd <- 1
        
        # loop through extrema which should be peaks and then valleys
        while(ii < len) {
            ii <- ii + 1 # this is a peak
            # reset peak finding if we had a peak and the next peak is bigger
            # than the last or the left min was small enough to reset.
            if(foundPeak && 
               (x[ii] > peakMag[length(peakMag)] || 
                leftMin < peakMag[length(peakMag)]-thresh)) {
                tempMag <- minMag
                foundPeak <- FALSE
            }
            
            #make sure we dont iterate past the length of our vector
            if(ii == len) {
                break
            }
            
            #found new peak that was larger than temp mag and threshold larger
            # than the minimum to its left
            if(x[ii] > tempMag && x[ii] > leftMin + thresh) {
                tempLoc <- ii
                tempMag <- x[ii]
            }
            
            ii <- ii + 1 # move onto the valley
            # come down at least thresh from peak
            if(!foundPeak && tempMag > thresh + x[ii]) {
                foundPeak <- TRUE
                leftMin <- x[ii]
                peakLoc[cInd] <- tempLoc
                peakMag[cInd] <- tempMag
                cInd <- cInd + 1
            } else if(x[ii] < leftMin) { #new left minima
                leftMin <- x[ii]
            }
        }
        
        # check end point
        if(x[length(x)] > tempMag &&
           x[length(x)] > leftMin + thresh) {
            peakLoc[cInd] <- len
            peakMag[cInd] <- x[length(x)]
            cInd <- cInd + 1
        } else if(!foundPeak && tempMag > minMag) { #check if we still need to add the last point
            peakLoc[cInd] <- tempLoc
            peakMag[cInd] <- tempMag
            cInd <- cInd + 1
        }
        
        # create output
        peakInds <- ind[peakLoc[1:(cInd-1)]]
        peakMag <- peakMag[1:(cInd-1)]
    } else {
        peakMag <- max(x)
        xInd <- which.max(x)
        if(peakMag > minMag + thresh) {
            peakInds <- inds[xInd]
        } else {
            peakMag <- NULL
            peakInds <- NULL
        }
    }
    
    # change sign of data if was finding minima
    if(extrema < 0) {
        peakMag <- -peakMgs
        x0 <- -x0
    }
    if(plot) {
        plot(1:length(x0), x0, type='l')
        points(peakInds, peakMag)
    }
    list(peakLoc=peakInds, peakMag=peakMag)
}

# called demod(c_half, mfreq, SampleRate, 'am')
# y signal, Fc carrier freq, Fs samplerate, method
#' @importFrom signal filtfilt butter
# PAMmisc
demod <- function(y, Fc, Fs, method='am') {
    if(Fc >= Fs/2) {
        stop('Invalid carrier frequency (higher than Nyquist)')
    }
    # TODO: apply to matrix if matrix here
    if(!is.null(dim(y))) {
        return()
    }
    
    len <- length(y) # change to be by matrix
    if(method %in% c('am', 'amdsb-sc', 'amdsb-tc', 'amssb')) {
        t <- seq(from=0, to=(len-1)/Fs, by=1/Fs)
        x_in <- y * cos(2*pi*Fc*t)
        bfilt <- butter(5, Fc*2/Fs)
        x_in <- filtfilt(bfilt, x_in)
        # x_in <- EndEffect(bfilt, x_in)
    }
    x_in
}

# EndEffect <- function(filt,x) {
#     signal::filtfilt(filt,c(rev(x),x,rev(x)))[(length(x) + 1):(2 * length(x))]
# }

# (var names copied from annamaria code)
# Za array depth
# R slant range from pamguard
# dt slant delay time
# c soundspeed
# SR perpendicular distance
# calculateDepth <- function(Za, R, dt, c, SR=NULL) {
# RENAME to SRE 
# PAMmisc
calculateDepth <- function(arrDepth, slantRange, delayTime, soundSpeed, perpDist=NULL) {
    extradist <- delayTime * soundSpeed
    # a <- slantRange + extradist
    # b <- slantRange
    # c <- arrDepth*2
    # 
    # rs <- b^2 + c^2
    # ls <- a^2 - rs
    # ls <- ls / (-2*b*c)
    # phi <- acos(ls)
    # 
    # depth <- slantRange*sin(phi-pi/2) + arrDepth
    # 
    # # min extra is 2 * arrDepth * minDepth/maxRange
    # depth2 <- ((slantRange+extradist)^2 - slantRange^2)/4/arrDepth
    extradist * (2*slantRange + extradist) / 4 / arrDepth
    # min/maxR * (slantRange + arrDepth*min/maxR) ~ minD * (slantRange/maxR)
    # depth
    # list(depth, depth2)
    # depth2
}

#' @importFrom tuneR readWave
#' @importFrom signal butter filter
# rename and in PAMmisc
doEchoPeak <- function(wav, filter=c(25e3, 70e3), sr=NULL, 
                       peakMin=.01, minTime=NULL, maxTime=NULL,
                       arrDepth=5.5, soundSpeed=1500,
                       maxRange=4e3, minDepth=500, 
                       plot=TRUE, channel=NULL, n=3,
                       clipLen=.03) {
    if(is.character(wav)) {
        wav <- readWave(wav, toWaveMC = TRUE, from=0, to=clipLen, units='seconds')
    }
    if(inherits(wav, 'WaveMC')) {
        sr <- wav@samp.rate
        if(is.null(channel)) {
            channel <- ncol(wav@.Data)
        }
        wav <- wav@.Data[, channel] / 2^(wav@bit - 1)
    }
    wav <- wav[1:(clipLen * sr)]
    lowFilt <- butter(8, filter[2]*2/sr, type='low')
    highFilt <- butter(8, filter[1]*2/sr, type='high')
    wav <- signal::filter(lowFilt, wav)
    wav <- signal::filter(highFilt, wav)
    # c <- acf(wav, demean=FALSE, lag.max=length(wav), plot=FALSE)$acf[, 1, 1]
    c <- ffwAcf(wav)
    mfreq <- mean(filter)
    z <- demod(c, mfreq, sr)
    peaks <- peakfinder(abs(z), peakMin)
    peaks$peakTime <- (peaks$peakLoc-1)/sr
    
    if(is.null(minTime)) {
        minTime <- 2*arrDepth*minDepth / (soundSpeed*maxRange)
    }
    if(is.null(maxTime)) {
        maxTime <- 2*arrDepth / soundSpeed
    }
    isWithin <- peaks$peakTime > minTime &
        peaks$peakTime < maxTime
    
    peaks$peakMag <- peaks$peakMag[isWithin]
    peaks$peakLoc <- peaks$peakLoc[isWithin]
    peaks$peakTime <- peaks$peakTime[isWithin]
    peakSort <- sort(peaks$peakMag, index.return=TRUE, decreasing=TRUE)
    
    topN <- peakSort$ix[1:n]
    if(plot) {
        # make sure to reset changed plotting settings
        oldPar <- par()
        on.exit(par(oldPar))
        if(all((par()$mfrow == 1))) {
            par(mfrow=c(2, 1))
            thisPlotIx <- c(1,1)
        } else {
            thisPlotIx <- par()$mfg[1:2]
        }
        # par(mfrow=c(2,1), mar=c(3, 2, 2, 1))
        par(mar=c(1.5, 2, 2, 1))
        plot(x=(0:(length(wav)-1))/sr, y=wav, type='l')
        wavPeak <- which.max(abs(wav))
        points(x=wavPeak/sr + peaks$peakTime[topN], y=rep(0, n), col='red', pch=16)
        points(x=wavPeak/sr + peaks$peakTime[topN[1]], y=c(0), col='green', pch=16)
        lines(x=rep(wavPeak/sr + maxTime, 2), y=c(-1, 1), col='blue', lty=3)
        
        par(mar=c(2, 2, 1.5, 1))
        par(mfg=thisPlotIx + c(1, 0))
        plot(x=0:(length(z)-1)/sr, y=abs(z), type='l', xlim=c(0, .02), xaxs='i', yaxs='i')
        points(x=peaks$peakTime[topN], y=peaks$peakMag[topN], col='red')
        points(x=peaks$peakTime[topN[1]], y=peaks$peakMag[topN[1]], col='green')
        lines(x=c(minTime, minTime), y=c(0, 1), lty=3, col='black')
        lines(x=c(maxTime, maxTime), y=c(0, 1), lty=3, col='blue')
        # if(!identical(par()$mfrow, c(2, 1))) {
        #     par(mfg=thisPlotIx + c(0, 1))
        # }
    }
    list(mag=peaks$peakMag[topN], time=peaks$peakTime[topN], wav=wav)
}

# processEchoWavFolder <- function(folder, nTop=5, 
#                                  peakMin=.01, arrDepth=5.5, filter=c(25e3, 70e3), 
#                                  minDepth=50, maxRange=6500,
#                                  progress=FALSE, nDo=Inf) {
#     files <- list.files(folder, pattern='wav$', full.names=TRUE, recursive=TRUE)
#     nDo <- min(nDo, length(files))
#     outTimes <- matrix(0, nrow=nDo, ncol=nTop)
#     outMags <- outTimes
#     if(progress) {
#         pb <- txtProgressBar(min=0, max=nDo, style=3)
#     }
#     for(i in 1:nDo) {
#         thisOut <- doEchoPeak(wav=files[i], peakMin = peakMin, arrDepth=arrDepth,
#                               minDepth=minDepth, maxRange=maxRange, 
#                               plot=FALSE, filter=filter, 
#                               n=nTop)
#         outTimes[i, ] <- thisOut$time
#         outMags[i, ] <- thisOut$mag
#         if(progress) {
#             setTxtProgressBar(pb, value=i)
#         }
#     }
#     list(times=outTimes, mags=outMags)
# }

#' @importFrom fftw FFT IFFT
# PAMmisc rename
ffwAcf <- function(x) {
    n <- length(x)
    fx <- FFT(c(x, rep(0, n)))
    x2 <- IFFT(abs(fx)^2)
    abs(x2[1:n]/x2[1])
}

# plotEchoWavs <- function(wavList, echo, filename, type=c('times', 'depth'), maxDepth=4e3,
#                          snrMin=5) {
#     # type <- match.arg(type)
#     n <- nrow(echo[[type[1]]])
#     evNames <- sapply(wavList[1:n], function(x) parseWavClipName(basename(x), 'Event'))
#     evIx <- split(1:n, evNames)
#     height <- 4
#     width <- 8 * length(type)
#     png(filename=filename, width=width, height=height * length(evIx), units='in', res=300)
#     on.exit(dev.off())
#     par(mfrow=c(length(evIx), length(type)))
#     for(i in seq_along(evIx)) {
#         for(j in seq_along(type)) {
#             yLims <- NULL
#             if(type[j] == 'depth') {
#                 yLims <- c(0, maxDepth)
#             }
#             if(type[j] == 'snr') {
#                 thisPlots <- echo[['times']][evIx[[i]], ]
#                 thisPlots[echo[['snr']][evIx[[i]]] < snrMin, ] <- NA
#                 plot(x=rep(seq_along(evIx[[i]]), ncol(echo[['times']])),
#                      y=thisPlots)
#                 points(x=seq_along(evIx[[i]]), y=thisPlots[, 1], col='red')
#                 # plot(x=seq_along(evIx[[i]]), y=echo[[type[j]]][evIx[[i]]], type='l')
#             } else {
#                 plot(x=rep(seq_along(evIx[[i]]), ncol(echo[[type[j]]])), 
#                      y=echo[[type[j]]][evIx[[i]], ], 
#                      main=names(evIx)[i], ylim=yLims)
#                 points(x=seq_along(evIx[[i]]), y=echo[[type[j]]][evIx[[i]], 1], col='red')
#             }
#         }
#     }
# }

# estSnr <- function(wav, len=.0025, sr=NULL, channel=NULL) {
#     if(is.character(wav)) {
#         wav <- readWave(wav, toWaveMC = TRUE)
#     }
#     if(inherits(wav, 'WaveMC')) {
#         sr <- wav@samp.rate
#         if(is.null(channel)) {
#             channel <- ncol(wav@.Data)
#         }
#         wav <- wav@.Data[, channel] / 2^(wav@bit - 1)
#     }
#     len <- round(len * sr)
#     signal <- PAMpal:::clipAroundPeak(wav, length=len, noise=FALSE)
#     noise <- PAMpal:::clipAroundPeak(wav, length=len, noise=TRUE)
#     10*log10(mean(signal^2)) - 10*log10(mean(noise^2))
# }

# PAMpal rename

calcDiveDepth <- function(study, wav, params, soundSpeed=1500, hpDepthError=1, clipLen=.03,
                          plot=TRUE, nPlot=400, nCol=5, outDir, locType='PGTargetMotion') {
    startTime <- Sys.time()
    if(length(wav) == 1 &&
       dir.exists(wav)) {
        wav <- list.files(wav, pattern='wav$', full.names=TRUE, recursive=TRUE)
    }
    if(isFALSE(plot)) {
        nPlot <- 0
    }
    wavExists <- file.exists(wav)
    if(!any(wavExists)) {
        stop('Could not find any of the wav files, check input')
    }
    if(any(!wavExists)) {
        warning(sum(!wavExists), ' wav files could not be found on disk.')
        wav <- wav[wavExists]
    }
    locData <- bind_rows(lapply(events(study), function(x) {
        localizations(x)[[locType]]
    }), .id='eventId')
    hasLoc <- locData$eventId[!is.na(locData$locLat)]
    # study <- study[hasLoc]
    clickData <- getClickData(study[hasLoc], measures=FALSE)
    if(!'hpDepth' %in% colnames(clickData)) {
        stop('Hydrophone depth data not found, please add with "addHyrophoneDepth"')
    }
    if(!all(c('Latitude', 'Longitude') %in% colnames(clickData))) {
        stop('GPS data not found, please add with "addGps"')
    }
    if(!dir.exists(outDir)) {
        dir.create(outDir)
    }
    ddCols <- c('UTC', 'eventId', 'Channel', 'locName', 'locLat', 'locLong',
                'perpDist', 'perpDistError', 'radialDist', 'hpDepth', 'Latitude',
                'Longitude', 'species')
    addedCols <- c('radialDist', 'wavFile', 
                   'maxTime', 'pair2Time', 'pair3Time',
                   'maxMag', 'pair2Mag', 'pair3Mag', 
                   'maxDepth', 'pair2Depth','pair3Depth')
    
    # dropping these if they already present for join
    clickData <- PAMpal:::dropCols(clickData, addedCols)
    
    clickData <- left_join(clickData, locData[c('eventId', 'locName', 
                                                'locLat', 'locLong',
                                                'perpDist', 'perpDistErr')],
                           by='eventId')
    clickData$radialDist <- distGeo(matrix(c(clickData$Longitude, clickData$Latitude), ncol=2),
                                    matrix(c(clickData$locLong, clickData$locLat), ncol=2))
    
    wavMatchDf <- data.frame(wavFile=wav,
                             UID=parseWavClipName(wav, 'UID'),
                             Channel=parseWavClipName(wav, 'channel'),
                             eventId=parseWavClipName(wav, 'event'))
    # drop this col in case already exists before join
    clickData$wavFile <- NULL
    clickData <- left_join(clickData, wavMatchDf, by=c('eventId', 'UID', 'Channel'))
    
    wavNoMatch <- !wav %in% clickData$wavFile[!is.na(clickData$wavFile)]
    if(sum(wavNoMatch) > 0) {
        warning(sum(wavNoMatch), ' wav files could not be matched to detections.')
        wav <- wav[!wavNoMatch]
    }
    hasMatch <- !is.na(clickData$wavFile) #clickData$wavFile %in% basename(wav)
    # clickData <- clickData[hasMatch, ]
    
    
    checkSpeciesParams(clickData$species, params)
    
    clickData <- split(clickData, clickData$eventId)
    clickData <- lapply(clickData, function(x) {
        x <- arrange(x, UTC)
        if(x$species[1] %in% names(params)) {
            thisParams <- params[[x$species[1]]]
        } else {
            thisParams <- params
        }
        # thisSr <- readWave(x$wavFile[1], header=TRUE)$sample.rate
        # lowFilt <- butter(8, thisParams$freqHigh*2/thisSr, type='low')
        # highFilt <- butter(8, thisParams$freqLow*2/thisSr, type='high')
        
        result <- vector('list', length=nrow(x))
        wavClips <- vector('list', length=nrow(x))
        nPlot <- min(nPlot, length(result))
        nRow <- ceiling(nPlot / nCol) * 2
        if(nPlot > 0) {
            png(file.path(outDir, paste0(x$eventId[1], '_Wavs.png')), width = 2*nCol, height=2*nRow, units='in', res=300)
            on.exit(dev.off())
            par(mfrow=c(nRow, nCol))
        }
        pb <- txtProgressBar(min=0, max=nrow(x), style=3)
        plotIx <- 1
        for(i in seq_along(result)) {
            if(is.na(x$wavFile[i])) {
                result[[i]] <- list(maxTime = NA,
                                    pair2Time = NA,
                                    pair3Time = NA,
                                    maxMag = NA,
                                    pair2Mag = NA,
                                    pair3Mag = NA,
                                    maxDepth = NA,
                                    pair2Depth = NA,
                                    pair3Depth = NA
                )
                next
            }
            doPlot <- plotIx <= nPlot
            if(doPlot) {
                plotCol <- ((plotIx-1) %% nCol) + 1
                plotRow <- floor((plotIx-1)/nCol)*2 + 1
                par(mfg=c(plotRow, plotCol))
            }
            thisWav <- readWave(x$wavFile[i], toWaveMC = TRUE, from=0, to=clipLen, units='seconds')
            #### POSSIBLE ADD DECIMATION ####
            # thisWav <- WaveMC(downsample(thisWav, 96e3), samp.rate=96e3, bit=16)
            waveHeightErr <- ifelse('waveHeight' %in% colnames(x), x$waveHeight[i], 0)
            thisEcho <- doEchoPeak(wav=thisWav,
                                   filter=c(thisParams$freqLow, thisParams$freqHigh),
                                   peakMin=.01,
                                   minTime=thisParams$minTime,
                                   maxTime=thisParams$maxTime,
                                   arrDepth=x$hpDepth[i] + hpDepthError + waveHeightErr,
                                   soundSpeed=soundSpeed,
                                   maxRange=thisParams$maxRange,
                                   minDepth=thisParams$minDepth,
                                   n=3,
                                   plot=doPlot,
                                   clipLen=clipLen)
            
            thisDepth <- calculateDepth(arrDepth=x$hpDepth[i], 
                                        slantRange =  x$radialDist[i],
                                        delayTime = thisEcho$time,
                                        soundSpeed = soundSpeed
            )
            outVals <- list(maxTime = thisEcho$time[1],
                            pair2Time = thisEcho$time[2],
                            pair3Time = thisEcho$time[3],
                            maxMag = thisEcho$mag[1],
                            pair2Mag = thisEcho$mag[2],
                            pair3Mag = thisEcho$mag[3],
                            maxDepth = thisDepth[1],
                            pair2Depth = thisDepth[2],
                            pair3Depth = thisDepth[3]
            )
            result[[i]] <- outVals
            wavClips[[i]] <- thisWav
            setTxtProgressBar(pb, value=i)
            plotIx <- plotIx + 1
        }
        result <- bind_rows(result)
        result <- cbind(x, result)
        if(nPlot > 0) {
            dev.off() # closing png plot for wav summary
            on.exit() # and removing that on.exit call
        }
        if(plot) {
            wavClips <- lapply(wavClips, function(x) {
                if(is.null(x)) return(x)
                #### HARD CODED ####
                x <- seewave::bwfilter(x, n=8, from=2e3, output='Wave')
                x <- downsample(x, samp.rate=9600)
                10*log10(abs(x@left)+1)
            })
            wavClips <- do.call(cbind, wavClips)
            png(file.path(outDir, paste0(x$eventId[1], '_Summary.png')), width=12, height=8, units='in', res=300)
            on.exit(dev.off())
            par(mfrow=c(2, 2))
            par(mar=c(2,4,2,1))
            # plot(x=seq_along(result$UTC), y=result$maxTime, col='red',
            #      main=paste0('Delay Times - ', x$eventId[1]))
            # points(x=seq_along(result$UTC), y=result$pair2Time, col='black')
            # points(x=seq_along(result$UTC), y=result$pair3Time, col='black')
            # result <- arrange(result, UTC)
            ## Depth Plot
            plot(x=result$UTC, y=-result$pair2Depth, col='black', ylim=c(-4e3,0),
                 main=paste0('Estimated Depth - ', x$eventId[1]), xlab='',
                 ylab='Depth (m)')
            points(x=result$UTC, y=-result$pair3Depth, col='black')
            points(x=result$UTC, y=-result$maxDepth, col='red')
            
            #### Echo-HARDCODED ####
            image(z=t(wavClips), x=1:ncol(wavClips), y=(1:nrow(wavClips))/9.600, 
                  xlab='', ylab='Time (ms)', main='Echogram',
                  col=viridisLite::viridis(32))
            par(mar=c(4, 4, 2, 1))
            ## Angle Plot
            plot(x=result$UTC, y=result$angle * 180 / pi,
                 main='Received Angle', ylab='Angle', xlab='UTC')
            
            ## Concat Clicks
            avgSpec <- calculateAverageSpectra(study[x$eventId[1]], evNum=1, plot=c(TRUE, FALSE))
        }
        result
    })
    
    clickData <- bind_rows(clickData)
    # dont need to carry these around
    locCols <- c('locName', 
                 'locLat', 'locLong',
                 'perpDist', 'perpDistErr')
    clickData <- PAMpal:::dropCols(clickData, locCols)
    endTime <- Sys.time()
    procTime <- round(as.numeric(difftime(endTime, startTime, units='secs')), 0)
    study <- detDataToStudy(study, clickData)
    cat('\nProcessing took ', procTime, ' seconds', sep='')
    study
    
}
# PAMpal already has this in writeEventClips, just adjust to work same as this
parseWavClipName <- function(x, mode=c('UID', 'UTC', 'event', 'channel')) {
    if(length(x) > 1) {
        return(sapply(x, function(i) parseWavClipName(i, mode), USE.NAMES = FALSE))
    }
    splitName <- strsplit(basename(x), '\\.')[[1]]
    switch(match.arg(mode),
           UID = {
               gsub('([0-9]*)CH[0-9_]*', '\\1', splitName[length(splitName)-1])
           },
           channel = {
               gsub('([0-9]*)CH([0-9])_.*', '\\2', splitName[length(splitName)-1])
           },
           UTC = {
               utcStr <- splitName[length(splitName)-1]
               milli <- as.numeric(substr(utcStr, stop=nchar(utcStr), start=nchar(utcStr)-2))/1e3
               utcStr <- substr(utcStr, stop=nchar(utcStr)-4, start=nchar(utcStr)-18)
               as.POSIXct(utcStr, format='%Y%m%d_%H%M%S', tz='UTC') + milli
           },
           event = {
               paste0(gsub('^Event_|^Detection_', '', splitName[1]),
                      '.',
                      splitName[2])
           },
           NA
    )
}

# matches the base wavfile name to the full name in "wav"
# PAMpal util to help above
matchWavFile <- function(df, wav) {
    if(length(unique(dirname(wav))) == 1) {
        tryWav <- paste0(dirname(wav[1]), '/', df$wavFile)
        if(all(file.exists(tryWav))) {
            df$wavFile <- tryWav
            return(df)
        }
    }
    matchDf <- data.frame(wavFile = basename(wav),
                          fullfile = wav)
    df <- left_join(df, matchDf, by='wavFile')
    if(sum(!is.na(df$fullfile)) == length(wav)) {
        df$wavFile <- df$fullfile
        df$fullfile <- NULL
        return(df)
    }
    for(i in 1:nrow(df)) {
        thisMatch <- wav[basename(wav) == df$wavFile[i]]
        if(length(thisMatch) != 1) {
            thisMatch <- NA
        }
        df$wavFile[i] <- thisMatch
    }
    df
}

# helper to check valid species parameter inputs
# PAMpal util to help above
checkSpeciesParams <- function(species, params) {
    reqCols <- c('freqLow', 'freqHigh', 'minDepth', 'maxRange')
    if(all(reqCols %in% names(params))) {
        return(TRUE)
    }
    species <- unique(species)
    if(any(is.na(species))) {
        stop('Species labels have not been assigned, either use "setSpecies" or',
             ' change input to list all of ', paste0(reqCols, collapse=', '))
    }
    matchSpecies <- species %in% names(params)
    if(!all(matchSpecies)) {
        noMatch <- species[!matchSpecies]
        warning('Species ', paste0(noMatch, collapse=', '), 
                ' were present in data but not in provided species params (',
                paste0(names(params), collapse=', '), ')')
    }
    speciesProper <- sapply(params, function(x) {
        all(reqCols %in% names(x))
    })
    if(!all(speciesProper)) {
        stop('Not all species have the required fields ', paste0(reqCols, collapse=', '))
    }
    TRUE
}

# filter by time/depth nearby to remove clicks that are impossible due
# to an animals max swim speed
# PAMpal rename to SRE/DD specific
filterClicks <- function(data, time=30, depth=NULL, speed=NULL,
                         maxDepth=4000, minCorr=.01) {
    if(is.null(depth) && is.null(speed)) {
        stop('Either depth or swim speed (m/s) must be provided')
    }
    if(is.null(depth)) {
        depth <- time * speed
    }
    if(is.AcousticStudy(data) || is.AcousticEvent(data)) {
        clicks <- getClickData(data)
    }
    if(is.character(data)) {
        if(!file.exists(data)) {
            stop('File ', data, ' does not exist')
        }
        clicks <- read.csv(data, stringsAsFactors = FALSE)
    }
    if(is.data.frame(data)) {
        clicks <- data
    }
    clicks$UTC <- PAMpal:::parseUTC(clicks$UTC)
    clicks <- split(clicks, clicks$eventId)
    clicks <- lapply(clicks, function(x) {
        filterBySpeed(x, time=time, depth=depth)
    })
    clicks <- bind_rows(clicks)
    clicks$keepClick[clicks[['maxDepth']] > maxDepth] <- FALSE
    clicks$keepClick[clicks[['maxMag']] < minCorr] <- FALSE
    if(is.AcousticStudy(data)) {
        data <- detDataToStudy(data, clicks)
        return(data)
    }
    clicks
}

# marks keepClick FALSE if moving too fast
# PAMpal util for above prob rename to swim speed?
filterBySpeed <- function(data, time=30, depth=NULL, speed=NULL) {
    if(is.null(depth) && is.null(speed)) {
        stop('Either depth or swim speed (m/s) must be provided')
    }
    if(is.null(depth)) {
        depth <- time * speed
    }
    # init values - keep only those non-NA
    if(!'keepClick' %in% colnames(data)) {
        data$keepClick <- !is.na(data$maxDepth)
    }
    if(all(is.na(data$maxDepth))) {
        return(data)
    }
    if(nrow(data) <= 1) {
        return(data)
    }
    for(i in seq_len(nrow(data))) {
        tDiffs <- abs(as.numeric(difftime(data$UTC[i], data$UTC[-i], units='secs')))
        dDiffs <- abs(data$maxDepth[i] - data$maxDepth[-i])
        # only keep if we wer alrady keeping and this test is still good
        data$keepClick[i] <- data$keepClick[i] & 
            isTRUE(any((tDiffs <= time) & (dDiffs <= depth)))
    }
    data
}

# PAMpal 
#' @import shiny
#' @import ggplot2
runDepthReview <- function(data) {
    inStudy <- is.AcousticStudy(data)
    inDf <- is.data.frame(data)
    if(!inStudy && !inDf) {
        stop('Input must be AcousticStudy or dataframe')
    }
    # using SHINYDATA as a temp var name to modify with <<- later from shiny
    if(inStudy) {
        SHINYDATA <- getClickData(data)
    }
    if(inDf) {
        SHINYDATA <- data
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
            SHINYDATA <- detDataToStudy(data, SHINYDATA)
        }
        return(invisible(SHINYDATA))
    })
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
                                label='Plot slant delay or depth',
                                choices=list('Slant delay'='maxTime', 'Depth'='maxDepth')),
                   actionButton('loadEcho', 
                                label='Load echogram data (slow)'),
                   actionButton("save","Save"), # when clicked saves brushed dataframe to file
                   actionButton('exit', 'Exit')
            ),
            column(width=9,
                   plotOutput(outputId = 'echogram')
            )
        )
    )
    
    # Server code begins here
    server <- function(input, output, session) {
        # making the dataset reactiveValues so that any changes in mt$data later could be reflected throughout
        smf <- reactiveValues(data=SHINYDATA,
                              echoData=list()) 
        # this populates dropdown with event names
        updateSelectInput(inputId='evSelect', choices=unique(SHINYDATA$eventId))
        yLabels <- list(
            'maxDepth' = 'Estimated depth (m)',
            'maxTime' = 'Measured slant delay (s)'
        )
        #### scatterplot ####
        output$scatterplot <- renderPlot({
            plotData <- smf$data %>% 
                dplyr::filter(eventId == input$evSelect)
            g <- ggplot(plotData, aes(x = UTC, y = .data[[input$plotValue]], col=keepClick)) +
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
            # if(all(is.na(plotData[[input$plotValue]]))) {
            #     g <- g +
            #         annotate('text', x=mean(plotData$UTC), y=1, label='All Values NA')
            # }
            # if(input$plotValue == 'maxDepth') {
            #     g <- g + ylim(-4000, 0)
            # }
            g
        })
        #### load echo data ####
        observeEvent(input$loadEcho, {
            # check if we've already loaded this events data
            if(input$evSelect %in% names(smf$echoData)) {
                thisEcho <- smf$echoData[[input$evSelect]]
            } else {
                # if we havent then load soundfiles
                showNotification('Loading wav clips...')
                thisEcho <- do.call(
                    cbind, 
                    lapply(smf$data$wavFile[smf$data$eventId == input$evSelect], function(x) {
                        if(is.na(x) || !file.exists(x)) {
                            return(NULL)
                        }
                        oneWav <- readWave(x, toWaveMC = TRUE, from=0, to=.03, units='seconds')
                        oneWav <- seewave::bwfilter(oneWav, n=8, from=2e3, output='Wave')
                        oneWav <- downsample(oneWav, samp.rate=9600)
                        10*log10(abs(oneWav@left)+1)
                    }))
                # store for easy access in future
                smf$echoData[[input$evSelect]] <- thisEcho
                showNotification('Done loading!')
            }
        })
        observeEvent(input$allFalse, {
            smf$data$keepClick[smf$data$eventId == input$evSelect] <- FALSE
        })
        #### echogram plot ####
        output$echogram <- renderPlot({
            if(input$evSelect %in% names(smf$echoData)) {
                toPlot <- smf$data$keepClick[smf$data$eventId == input$evSelect & !is.na(smf$data$wavFile)]
                if(any(toPlot)) {
                    echoData <- smf$echoData[[input$evSelect]]
                    echoData <- echoData[, toPlot]
                    par(mar=c(3.1, 4.1, 3.1, 2.1))
                    image(z=t(echoData), x=1:ncol(echoData), y=(1:nrow(echoData))/9.600, 
                          xlab='', ylab='Time (ms)', main='Echogram',
                          col=viridis32)
                }
            }
        })
        #### observe brush ####
        observe({
            df = brushedPoints(smf$data, brush = input$plot_brush_, yvar=input$plotValue, allRows = TRUE) # get column with false and true with AllRows = T
            # removing datapoints selected by brush, select values are stored as chars
            smf$data$keepClick[df$selected_ & df$eventId == input$evSelect] <- input$paintFlag == 'TRUE' 
        })
        #### save progress ####
        observeEvent(input$save, {
            # <<- puts it back to environment of original function call
            SHINYDATA <<- smf$data
            # write.csv(smf$data, paste(SaveFolder,Subset_metadata$eventId[1],"_filtered_Click_metadata_&_dts.csv",sep=""), row.names = FALSE) 
            showNotification('Saved!')
        })
        observeEvent(input$exit, {
            stopApp()
        })
    }
    
    # Create a Shiny app object
    app <- shinyApp(ui = ui, server = server)
    
    #open shiny app
    runApp(app)  
}
# output from viridisLite::viridis(32)
viridis32 <- c("#440154FF", "#470D60FF", "#48196BFF", "#482475FF", "#472E7CFF",
               "#453882FF", "#424186FF", "#3E4B8AFF", "#3A548CFF", "#365D8DFF",
               "#32658EFF", "#2E6D8EFF", "#2B758EFF", "#287D8EFF", "#25848EFF",
               "#228C8DFF", "#1F948CFF", "#1E9C89FF", "#20A386FF", "#25AB82FF",
               "#2EB37CFF", "#3ABA76FF", "#48C16EFF", "#58C765FF", "#6ACD5BFF",
               "#7ED34FFF", "#93D741FF", "#A8DB34FF", "#BEDF26FF", "#D4E21AFF",
               "#E9E51AFF", "#FDE725FF")

# PAMpal 
addWaveHeight <- function(x, height) {
    if(is.AcousticStudy(x)) {
        for(e in seq_along(events(x))) {
            for(d in seq_along(detectors(x[[e]]))) {
                thisType <- attr(x[[e]][[d]], 'calltype')
                x[[e]][[d]] <- addWaveHeight(x[[e]][[d]], height)
                attr(x[[e]][[d]], 'calltype') <- thisType
            }
        }
        return(x)
    }
    if(is.data.frame(x)) {
        if(is.numeric(height)) {
            x$waveHeight <- height
            return(x)
        }
        if(!('UTC' %in% colnames(height))) {
            stop('Wave height data must have column "UTC"')
        }
        if(!inherits(height$UTC, 'POSIXct')) {
            height$UTC <- PAMpal:::parseUTC(height$UTC)
        }
        if(sum(c('waveHeight', 'beaufort') %in% colnames(height)) < 1) {
            stop('Wave height data must have column "waveHeight" or "beaufort"')
        }
        setDT(height)
        if(!('waveHeight' %in% colnames(height))) {
            setkeyv(height, 'beaufort')
            height <- PAMpal:::ppVars()$bftHeight[height, roll=-Inf]
        }
        setkeyv(height, 'UTC')
        setDT(x)
        setkeyv(x, 'UTC')
        x <- height[x, roll=12*60*60]
        setDF(x)
        return(x)
    }
}

#' @importFrom readxl read_excel
# get effort, match bft to wave height, then set key for matching later
# this should probably not be a built in, have people load their own shit
# as a df and do renaming
readBftNEFSC <- function(x) {
    if(is.character(x)) {
        if(grepl('\\.xls', x)) {
            effort <- read_excel(x)
            effort <- select(effort, c('BEAUFORT', 'BEGDATETIMEGMT'))
        } else if(grepl('\\.csv', x)) {
            effort <- read.csv(x, stringsAsFactors = FALSE)
            effort <- select(effort, c('BEAUFORT', 'BEGDATETIMEGMT'))
        }
    }
    if(is.character(effort[['BEGDATETIMEGMT']])) {
        effort[['BEGDATETIMEGMT']] <- PAMpal:::parseUTC(effort[['BEGDATETIMEGMT']],
                                                        format = c('%d-%b-%y %H:%M:%S',
                                                                   '%m/%d/%Y %H:%M:%S',
                                                                   '%Y-%m-%d %H:%M:%S'))
    }
    effort <- rename(effort, 'UTC' = 'BEGDATETIMEGMT', 'beaufort'='BEAUFORT')
    effort
}

# goes from a getXXXData det dataframe back into the study it came from
# already pal'd
detDataToStudy <- function(study, dets) {
    colsToDrop <- c('eventId', 'detectorName', 'db', 'species')
    measNames <- names(getMeasures(study))
    colsToDrop <- unique(c(colsToDrop, measNames))
    dets <- split(dets, dets[['eventId']])
    detEvs <- names(dets)
    for(e in seq_along(dets)) {
        if(!names(dets)[e] %in% names(events(study))) {
            next
        }
        thisDet <- split(dets[[e]], dets[[e]][['detectorName']])
        thisDet <- lapply(thisDet, function(x) {
            PAMpal:::dropCols(x, colsToDrop)
        })
        for(d in seq_along(thisDet)) {
            ct <- attr(study[[names(dets)[e]]][[names(thisDet)[d]]], 'calltype')
            attr(thisDet[[d]], 'calltype') <- ct
            study[[names(dets)[e]]][[names(thisDet)[d]]] <- thisDet[[d]]
        }
    }
    study
}

# dt error, array depth error, slant range error (all pcts, then sum sqr