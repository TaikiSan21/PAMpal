#' @title Calculate Depth from Echoes
#'
#' @description Calculate the estimated depth of echolocation clicks using
#'   surface reflected echoes. This uses the time delay between the received
#'   signal and its surface echo to estimate the depth of a calling animal.
#'   Requires that a set of waveform clips has been
#'   created using \link{writeEventClips}, and that events have been localized.
#'
#' @param x \linkS4class{AcousticStudy} object
#' @param wav either folder containing wave clips or list of wave clip files
#' @param clipLen length (seconds) of clip to analyze
#' @param spParams list of species-specific parameters, see details
#' @param soundSpeed sound speed (meters/second) to use for calculations
#' @param hpDepthError maximum error (meters) in hydrophone depth measurement
#' @param locType name of localization, note that this function is not computing
#'   any localization, only using previously calculated
#' @param plot logical flag to create summary plots
#' @param nPlot number of waveform plots to create for summary
#' @param nCol number of columns for waveform summary plot
#' @param plotDir directory to store plot outputs, default \code{NULL} will result
#'   in no plots being created
#' @param progress logical flag to show progress bar
#' @param verbose logical flag to show messages
#'
#' @author Taiki Sakai \email{taiki.sakai@@noaa.gov}
#'
#' @details \code{spParams} allows for species-specific filtering and acceptable
#'   echo time delays to be specified. These are provided as a list with elements
#'   \code{freqLow} and \code{freqHigh} specifying the lower and upper ends of a
#'   bandpass filter to apply to the signals (in Hz), which can aid in properly
#'   detecting the echoes. Parameters \code{minTime} and \code{maxTime} can
#'   also be supplied to define ranges on allowed time delay values. Alternatively
#'   if \code{maxTime} is \code{NULL} or not present it will be calculated from the
#'   hydrophone geometry, and \code{minTime} can be calculated from geometry by
#'   providing \code{minDepth} and \code{maxRange} as the minimum detectable depth
#'   and maximum detectable range (in meters).
#'
#'   If the same values for these parameters should be used for all detections in
#'   \code{x}, then \code{spParams} can be provided as a list with each parameter named,
#'   e.g. \preformatted{
#'   list(freqLow=10e3, freqHigh=50e3, minTime=.001, maxTime=NULL)}
#'
#'   If different values should be used for different species, then \code{spParams}
#'   must be a named list where the names match the species in \code{x}, providing
#'   a separate list of values for each species. e.g.
#'   \preformatted{
#'   list(Zc=list(freqLow=10e3, freqHigh=50e3, minTime=.001, maxTime=NULL),
#'        Pm=list(freqLow=2e3, freqHigh=16e3, minTime=.001, maxTime=NULL))}
#'
#' @return the AcousticStudy object \code{x} with estimated dive depth outputs
#'   added for each detection that had a matching wav clip file in \code{wav}.
#'   Detections that either did not have matchinf wav files or did not have
#'   localizations will have \code{NA} for all dive depth outputs. The depth outputs
#'   are
#'   \describe{
#'      \item{maxTime}{Delay time with maximum correlation value}
#'      \item{pair2Time}{Delay time with second highest correlation value}
#'      \item{pair3Time}{Delay time with third highest correlation value}
#'      \item{maxMag}{Correlation magnitude for "maxTime"}
#'      \item{pair2Mag}{Correlation magnitude for "pair2Time"}
#'      \item{pair3Mag}{Correlation magnitude for "pair3Time"}
#'      \item{maxDepth}{Calculated depth for "maxTime"}
#'      \item{pair2Depth}{Calculated depth for "pair2Time"}
#'      \item{pair3Depth}{Calculated depth for "pair3Time"}
#'   }
#'
#' @examples
#' # example not run because it requires access to large files not present
#' # in the package testing material
#' \dontrun{
#' study <- addRecordings(study, folder='path/to/recordings')
#' wavPath <- 'path/to/wavFiles'
#' writeEventClips(study, outDir=wavPath, mode='detection')
#' study <- calculateEchoDepth(study, wav=wavPath)
#' }
#'
#' @importFrom PAMmisc findEchoTimes
#' @importFrom graphics layout
#' @importFrom grDevices dev.off png
#'
#' @export
#'
calculateEchoDepth <- function(x,
                               wav,
                               clipLen=.03,
                               spParams=NULL,
                               soundSpeed=1500,
                               hpDepthError=1,
                               locType='PGTargetMotion',
                               plot=TRUE,
                               nPlot=400,
                               nCol=5,
                               plotDir=NULL,
                               progress=TRUE,
                               verbose=TRUE) {
    startTime <- Sys.time()
    if(length(wav) == 1 &&
       dir.exists(wav)) {
        wav <- list.files(wav, pattern='wav$', full.names=TRUE, recursive=TRUE)
    }
    if(isFALSE(plot)) {
        nPlot <- 0
        plotDir <- NULL
    }
    if(isTRUE(plot) && is.null(plotDir)) {
        stop('Must specify output directory "plotDir" in order to create summary plots.')
    }
    wavExists <- file.exists(wav)
    if(!any(wavExists)) {
        stop('Could not find any of the wav files, check input')
    }
    if(any(!wavExists)) {
        warning(sum(!wavExists), ' wav files could not be found on disk.')
        wav <- wav[wavExists]
    }
    locData <- bind_rows(lapply(events(x), function(x) {
        localizations(x)[[locType]]
    }), .id='eventId')
    hasLoc <- locData$eventId[!is.na(locData$locLat)]
    # x <- x[hasLoc]
    clickData <- getClickData(x[hasLoc], measures=FALSE)
    if(!'hpDepth' %in% colnames(clickData)) {
        stop('Hydrophone depth data not found, please add with "addHydrophoneDepth"')
    }
    if(!all(c('Latitude', 'Longitude') %in% colnames(clickData))) {
        stop('GPS data not found, please add with "addGps"')
    }

    addedCols <- c('radialDist', 'wavFile',
                   'maxTime', 'pair2Time', 'pair3Time',
                   'maxMag', 'pair2Mag', 'pair3Mag',
                   'maxDepth', 'pair2Depth','pair3Depth')

    # dropping these if they already present for join
    clickData <- dropCols(clickData, addedCols)

    clickData <- left_join(clickData, locData[c('eventId', 'locName',
                                                'locLat', 'locLong',
                                                'perpDist', 'perpDistErr')],
                           by='eventId')
    clickData$radialDist <- distGeo(matrix(c(clickData$Longitude, clickData$Latitude), ncol=2),
                                    matrix(c(clickData$locLong, clickData$locLat), ncol=2))

    wavMatchDf <- data.frame(wavFile=wav,
                             UID=parseEventClipName(wav, 'UID'),
                             Channel=parseEventClipName(wav, 'channel'),
                             eventId=parseEventClipName(wav, 'event'))
    # drop this col in case already exists before join
    clickData$wavFile <- NULL
    clickData <- left_join(clickData, wavMatchDf, by=c('eventId', 'UID', 'Channel'))

    wavNoMatch <- !wav %in% clickData$wavFile[!is.na(clickData$wavFile)]
    hasMatch <- !is.na(clickData$wavFile) #clickData$wavFile %in% basename(wav)
    if(!any(hasMatch)) {
        stop('No detections had matching wav files, check directory.')
    }
    if(sum(wavNoMatch) > 0) {
        warning(sum(wavNoMatch), ' wav files could not be matched to detections.')
        wav <- wav[!wavNoMatch]
    }
    # clickData <- clickData[hasMatch, ]

    spParams <- checkSpeciesParams(clickData$species, spParams)

    clickData <- split(clickData, clickData$eventId)
    if(progress) {
        if(verbose) {
            cat('Processing', length(clickData), 'events...')
        }
        evIx <- 1
    }
    if(!is.null(plotDir) &&
       !dir.exists(plotDir)) {
        dir.create(plotDir)
    }
    clickData <- lapply(clickData, function(ev) {
        ev <- arrange(ev, UTC)
        if(ev$species[1] %in% names(spParams)) {
            thisParams <- spParams[[ev$species[1]]]
        } else {
            thisParams <- spParams
        }
        result <- vector('list', length=nrow(ev))
        wavClips <- vector('list', length=nrow(ev))
        hasWav <- !is.na(ev$wavFile)
        nPlot <- min(nPlot, sum(hasWav))
        nRow <- ceiling(nPlot / nCol) * 2
        if(nPlot > 0) {
            png(file.path(plotDir, paste0(ev$eventId[1], '_Wavs.png')),
                width = 2*nCol, height=2*nRow, units='in', res=300)
            on.exit(dev.off())
            par(mfrow=c(nRow, nCol))
        }
        if(progress) {
            cat(paste0('\n', evIx, '/', length(clickData), '\n'))
            pb <- txtProgressBar(min=0, max=nrow(ev), style=3)
            evIx <<- evIx + 1
        }
        plotIx <- which(hasWav)[round(seq(from=1, to=sum(hasWav), length.out=nPlot), 0)]
        for(i in seq_along(result)) {
            if(is.na(ev$wavFile[i])) {
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
            doPlot <- i %in% plotIx
            if(doPlot) {
                plotCol <- ((which(i==plotIx)-1) %% nCol) + 1
                plotRow <- floor((which(i==plotIx)-1)/nCol)*2 + 1
                par(mfg=c(plotRow, plotCol))
            }
            #### POSSIBLE ADD DECIMATION ####
            # thisWav <- WaveMC(downsample(thisWav, 96e3), samp.rate=96e3, bit=16)
            waveHeightErr <- ifelse('waveHeight' %in% colnames(ev), ev$waveHeight[i], 0)

            if(is.null(thisParams$minTime)) {
                minTime <- 2*(ev$hpDepth[i] + hpDepthError + waveHeightErr) * thisParams$minDepth /
                    (soundSpeed * thisParams$maxRange)
            } else {
                minTime <- thisParams$minTime
            }
            if(is.null(thisParams$maxTime)) {
                maxTime <- 2 * (ev$hpDepth[i] + hpDepthError + waveHeightErr) / soundSpeed
            } else {
                maxTime <- thisParams$maxTime
            }
            thisEcho <- findEchoTimes(wav=ev$wavFile[i],
                                      filter=c(thisParams$freqLow, thisParams$freqHigh),
                                      peakMin=.01,
                                      minTime=minTime,
                                      maxTime=maxTime,
                                      n=3,
                                      plot=doPlot,
                                      clipLen=clipLen,
                                      plotText=paste0('UID: ',ev$UID[i], '\nIndex: ', which(i==which(hasWav))))

            thisDepth <- calculateDepth(arrDepth=ev$hpDepth[i],
                                        slantRange =  ev$radialDist[i],
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
            wavClips[[i]] <- thisEcho$wav
            if(progress) {
                setTxtProgressBar(pb, value=i)
            }
        }
        result <- bind_rows(result)
        result <- cbind(ev, result)
        if(nPlot > 0 &&
           !is.null(plotDir)) {
            dev.off() # closing png plot for wav summary
            on.exit() # and removing that on.exit call
        }
        if(plot) {
            png(file.path(plotDir, paste0(ev$eventId[1], '_Summary.png')), width=12, height=8, units='in', res=300)
            on.exit(dev.off())
            f <- layout(
                matrix(c(1,1,2,2,3,3, 4,4,4,5,5,5), ncol=2, byrow=F)#,
            )
            par(mar=c(2,4,2,1))
            result <- arrange(result, UTC)

            plot(x=(result$UTC), y=result$pair2Time * 1e3, col='black',
                 main=ev$eventId[1],
                 ylab='Time (ms)')
            points(x=(result$UTC), y=result$pair3Time * 1e3, col='black')
            points(x=(result$UTC), y=result$maxTime * 1e3, col='red')
            ## Depth Plot
            plot(x=result$UTC, y=-result$pair2Depth, col='black', ylim=c(-4e3,0),
                 main='^ Delay Times ^ - v Estimated Depth v',
                 xlab='', ylab='Depth (m)')
            points(x=result$UTC, y=-result$pair3Depth, col='black')
            points(x=result$UTC, y=-result$maxDepth, col='red')

            ## Angle Plot
            plot(x=result$UTC, y=result$angle * 180 / pi,
                 main='Received Angle', ylab='Angle', xlab='UTC')

            ## echogram
            thisSr <- readWave(wav[1], header=TRUE)$sample.rate
            wavClips <- do.call(cbind, wavClips)
            plotEchogram(wavClips, q=c(.01, .999), sr=thisSr)
            par(mar=c(4, 4, 2, 1))

            ## Concat Clicks
            avgSpec <- calculateAverageSpectra(x[ev$eventId[1]], evNum=1, plot=c(TRUE, FALSE))
        }
        result
    })

    clickData <- bind_rows(clickData)
    # dont need to carry these around
    locCols <- c('locName',
                 'locLat', 'locLong',
                 'perpDist', 'perpDistErr')
    clickData <- dropCols(clickData, locCols)
    endTime <- Sys.time()
    procTime <- round(as.numeric(difftime(endTime, startTime, units='secs')), 0)
    x <- detDataToStudy(x, clickData)
    if(verbose) {
        cat('\nProcessing took ', procTime, ' seconds', sep='')
    }
    x

}

calculateDepth <- function(arrDepth, slantRange, delayTime, soundSpeed, perpDist=NULL) {
    # a <- slantRange + extradist
    # b <- slantRange
    # c <- arrDepth*2
    # rs <- b^2 + c^2
    # ls <- a^2 - rs
    # ls <- ls / (-2*b*c)
    # phi <- acos(ls)
    # depth <- slantRange*sin(phi-pi/2) + arrDepth
    extradist <- delayTime * soundSpeed
    extradist * (2*slantRange + extradist) / 4 / arrDepth
}

plotEchogram <- function(wavMat, sr, q=c(.01, .999)) {
    wavMat <- 10*log10(abs(wavMat) + 1e-10)
    lims <- quantile(wavMat, q)
    wavMat[wavMat < lims[1]] <- lims[1]
    wavMat[wavMat > lims[2]] <- lims[2]
    viridis32 <- c("#440154FF", "#470D60FF", "#48196BFF", "#482475FF", "#472E7CFF",
                   "#453882FF", "#424186FF", "#3E4B8AFF", "#3A548CFF", "#365D8DFF",
                   "#32658EFF", "#2E6D8EFF", "#2B758EFF", "#287D8EFF", "#25848EFF",
                   "#228C8DFF", "#1F948CFF", "#1E9C89FF", "#20A386FF", "#25AB82FF",
                   "#2EB37CFF", "#3ABA76FF", "#48C16EFF", "#58C765FF", "#6ACD5BFF",
                   "#7ED34FFF", "#93D741FF", "#A8DB34FF", "#BEDF26FF", "#D4E21AFF",
                   "#E9E51AFF", "#FDE725FF")
    image(z=t(wavMat), x=1:ncol(wavMat), y=(1:nrow(wavMat))/sr*1e3,
          xlab='', ylab='Time (ms)', main='Echogram',
          col=viridis32)
}

checkSpeciesParams <- function(species, params) {
    if(is.null(params)) {
        warning('No "spParams" specified. No filtering will be applied to',
                ' signals, and allowed delay times may be inaccurate.')
        return(list(freqLow = 0, freqHigh=NULL, minTime=.001, maxTime=NULL))
    }
    if(hasParamNames(params)) {
        return(params)
    }

    species <- unique(species)
    if(any(is.na(species))) {
        stop('Species labels have not been assigned, either use "setSpecies" or',
             ' change "spParams" to only list required parameters.')
    }
    matchSpecies <- species %in% names(params)
    if(!all(matchSpecies)) {
        noMatch <- species[!matchSpecies]
        warning('Species ', paste0(noMatch, collapse=', '),
                ' were present in data but not in provided species params (',
                paste0(names(params), collapse=', '), ')')
    }
    speciesProper <- sapply(params, function(x) {
       hasParamNames(x)
    })
    if(!all(speciesProper)) {
        stop('Not all species have the required fields, check "spParams".')
    }
    params
}

hasParamNames <- function(x) {
    all(c('freqLow', 'freqHigh') %in% names(x)) &&
        (all(c('minTime', 'maxTime') %in% names(x)) ||
             all(c('minDepth', 'maxRange') %in% names(x)))
}
