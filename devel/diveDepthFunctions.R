library(geosphere)
library(tuneR)
library(data.table)
library(PAMpal)
library(readxl)
# Updated 06-03-2021 TNS
# 5-20 changed to ArrayDepthPercentErr
# 6-30 Updated plotting code to show event name and also groups PNGs by event
# instead of lumping all plots into one file
# 12-15-2021 Updated plotting code to store in list instead of matrix so that
# 12-20-2021 Updated to skip clips too short
# different length wav clips do not cause crash

export_diveDepthNEFSC <- function(x, outDir=NULL, file=NULL, wavFolder, waveHeight,
                                  filter=NULL, clipLength=.03, soundSpeed=1500, depthSensAcc=1,
                                  plot=TRUE, progress=TRUE, overwrite=TRUE) {
    hasLocalization <- sapply(events(x), function(e) {
        tarMo <- localizations(e)$PGTargetMotion
        if(is.null(tarMo)) {
            return(FALSE)
        }
        !is.na(tarMo$TMLatitude1) &
            !is.na(tarMo$TMPerpendicularDistance1)
    })
    if(!any(hasLocalization)) {
        stop('No events had Target Motion Localizations')
    }
    if(any(!hasLocalization)) {
        warning(sum(!hasLocalization), ' out of ', length(hasLocalization),
                ' events did not have Target Motion Localizations')
    }
    x <- x[hasLocalization]
    if(!dir.exists(wavFolder)) {
        dir.create(wavFolder)
    }
    wavFiles <- list.files(wavFolder, pattern='\\.wav$', full.names = TRUE)
    if(overwrite ||
       length(wavFiles) == 0) {
        if(progress) {
            cat('Writing wav clips...\n')
            pb <- txtProgressBar(min=0, max=length(events(x)), style=3)
        }
        for(e in seq_along(events(x))) {
            thisFilt <- parseFilter(x[[e]], filter)
            wavFiles <- c(wavFiles, writeEventClips(x[e],
                                                    buffer = c(0, clipLength),
                                                    outDir = wavFolder,
                                                    mode='detection',
                                                    filter = thisFilt,
                                                    useSample = TRUE,
                                                    progress = FALSE))
            if(progress) {
                setTxtProgressBar(pb, value=e)
            }
        }
    }
    ddTable <- createDiveDepthTable(x, waveHeight=waveHeight, soundSpeed=soundSpeed,
                                    depthSensAcc=depthSensAcc, wavFolder=wavFolder,
                                    clipLength = clipLength)
    if(is.null(outDir)) {
        outDir <- wavFolder
    }

    clipList <- plotEchoSummary(wavFolder, outDir=outDir, ncol=5, nrow=100,
                                start=1, end=clipLength,
                                species=x, arrayDepth=ddTable, soundSpeed=soundSpeed)

    if(is.null(file)) {
        file <- paste0('DiveDepthMetadata_', id(x), '.csv')
    }
    write.csv(ddTable, row.names=FALSE, file=file.path(outDir, file))
    list(wavList = clipList, table = ddTable, file=file.path(outDir, file))
}
# x is AcSt
parseFilter <- function(x, filter) {
    if(is.null(filter)) {
        return(0)
    }
    if(is.numeric(filter)) {
        return(filter)
    }
    if(is.data.frame(filter)) {
        sp <- gsub('[[:punct:]]', '', tolower(species(x)$id))
        if(sp %in% filter$species) {
            from <- filter$from[filter$species == sp]
            if(is.na(from)) {
                from <- 0
            }
            to <- filter$to[filter$species == sp]
            if(is.na(to)) {
                return(from)
            }
            return(c(from, to))
        }
    }
    NULL
}

#' @importFrom readxl read_excel
# get effort, match bft to wave height, then set key for matching later
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

addWaveHeight <- function(x, height) {
    if(is.numeric(height)) {
        x$waveHeight <- height
        return(x)
    }
    if(!('UTC' %in% colnames(height))) {
        stop('Wave height data must have column "UTC"')
    }
    if(!inherits(height$UTC, 'POSIXct')) {
        stop('"UTC" must be converted to POSIXct')
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
    setkeyv(x, 'ClickTime')
    x <- height[x, roll=12*60*60]
    setDF(x)
    x
}

####

# algs
# least squares
# 2d simplex optimisation
# 3d simplex optimisiation
# mcmc localisation
####
createDiveDepthTable <- function(x, waveHeight, soundSpeed=1500, depthSensAcc = 1, wavFolder, clipLength) {
    clickData <- getClickData(x)
    # browser()
    if(!all(c('Longitude', 'Latitude') %in% colnames(clickData))) {
        stop('Lat/Long data has not been added, use "addGps" first.')
    }
    if(!('hpDepth' %in% colnames(clickData))) {
        stop('Hydrophone depth data has not been added, use "addHydrophoneDepth" first.')
    }
    clickData <- select(clickData, c('UID', 'UTC', 'species', 'Latitude', 'Longitude', 'hpDepth', 'eventId')) %>%
        distinct()
    tmData <- bind_rows(lapply(events(x), function(e) {
        c(id = id(e),
          localizations(e)$PGTargetMotion[c('locLat', 'locLong',
                                            'perpDist', 'perpDistErr')])
    }))
    clickData <- left_join(clickData, tmData, by = c('eventId' = 'id'))
    clickData$RadialDistM <- distGeo(matrix(c(clickData$Longitude, clickData$Latitude), ncol=2),
                                     matrix(c(clickData$locLong, clickData$locLat), ncol=2))
    clickData <- select(clickData,
                        UID,
                        EventId = eventId,
                        GroupLat = locLat,
                        GroupLon = locLong,
                        ShipLat = Latitude,
                        ShipLon = Longitude,
                        SlantRange = perpDist,
                        SlantRangeErr = perpDistErr,
                        Species = species,
                        ClickTime = UTC,
                        ArrayDepth = hpDepth,
                        RadialDistM) %>%
        mutate(SlantRangePercentErr = SlantRangeErr/SlantRange * 100)
    # this currently matches each detection time to the previous effort start time
    clickData <- addWaveHeight(clickData, waveHeight)
    clickData <- rename(clickData, 'Waveheight' = 'waveHeight')
    clickData$DepthSensAccuracy <- depthSensAcc
    clickData$SoundSpeed <- soundSpeed
    clickData$ArrayDepthPercentErr <- (clickData$Waveheight + clickData$DepthSensAccuracy) / clickData$ArrayDepth * 100
    if(length(wavFolder) == 1 &&
       dir.exists(wavFolder)) {
        wavFiles <- list.files(wavFolder, pattern= '\\.wav', full.names = TRUE)
    } else if(any(file.exists(wavFolder))) {
        wavFiles <- wavFolder
    } else {
        stop('No wav files found')
    }

    wavLen <- unname(sapply(wavFiles, function(w) {
        hdr <- readWave(w, header = TRUE)
        hdr$samples / hdr$sample.rate
    }))
    wavDf <- data.frame(Filename = basename(wavFiles),
                        UID = parseWavClipName(basename(wavFiles), 'UID'),
                        wavLen)
    clickData <- left_join(wavDf, clickData, by = c('UID' = 'UID'))

    tooShort <- clickData$wavLen < .99 * clipLength
    if(any(tooShort)) {
        warning('Removing ', sum(tooShort), ' wav clips that were shorter than the ',
                'desired clip length.', call.=FALSE)
        clickData <- clickData[!tooShort, ]
    }
    select(clickData,
           Filename,
           EventId,
           UID,
           ClickTime,
           Species,
           SoundSpeed,
           # Beaufort,
           Waveheight,
           ArrayDepth,
           DepthSensAccuracy,
           ArrayDepthPercentErr,
           GroupLat,
           GroupLon,
           ShipLat,
           ShipLon,
           RadialDistM,
           SlantRange,
           SlantRangeErr,
           SlantRangePercentErr)
}

##### PARAMS ####
# dir - folder of wav clips created by writeEventClips
# outDir - folder to create plot in
# start - starting point of clip to draw, in samples
# end - ending point of clip to draw, in samples. If NULL will use entire wav file
# from - lower limit for bandpass filter, in Hz
# to - upper end for bandpass filter, in Hz
# nrow, ncol - dimensions for grid of wav clips to create. If only one is specified,
#    other dimension will be determined to fit all clips in one large image. If both
#    are specified, then multiple PNGs will be created  with those dimensions
# vlines - locations to draw dashed vertical lines, in seconds
# species - One of Cuviers, Trues, Gervais, MmMe, Sowerbys, Blainvilles. Any punctation
#    will be removed and all characters set to lowercase before matching. These are used to
#    set specific filter levels.
# soundSpeed - default 1500 m/s
# arrayDepth - depth of array in meters. Used in combination with soundSpeed
#    to determine maximum echo time. Can also be an AcousticStudy object, in which
#    case matching depth data will be retrieved for each individual click
plotEchoSummary <- function(dir, outDir='.',
                            start=1, end=NULL,
                            from=0, to=NULL,
                            nrow=NULL, ncol=NULL,
                            vlines=NULL,
                            species=NULL,
                            soundSpeed=1500, arrayDepth = 15) {
    if(!dir.exists(dir)) {
        stop('Could not locate folder ', dir)
    }
    wavFiles <- list.files(dir, pattern='\\.wav$', full.names = TRUE)

    if(is.null(end)) {
        end <- readWave(wavFiles[1], header=TRUE)$samples
    }
    # check if in seconds
    if(is.numeric(end) &&
       end < 10) {
        end <- end * readWave(wavFiles[1], header=TRUE)$sample.rate
    }
    # dbMat <- matrix(0, nrow=length(signal::decimate(start:end, factor)), ncol=length(wavFiles))

    # wavMat <- matrix(0, nrow=(length(start:end)), ncol=length(wavFiles))
    wavList <- vector('list', length=length(wavFiles))

    cat('Reading wav files...\n')
    pb <- txtProgressBar(min=0, max=length(wavFiles), style=3)
    tooShort <- numeric(0)
    for(w in seq_along(wavFiles)) {
        setTxtProgressBar(pb, value = w)
        thisWav <- readWave(wavFiles[w], from=start, to=end)
        thisWav@left <- thisWav@left - mean(thisWav@left)
        thisWav <- thisWav / 2^(thisWav@bit-1)
        thisSr <- thisWav@samp.rate
        if(length(thisWav@left) < .99*(end - start + 1)) {
            tooShort <- c(tooShort, w)
            next
        }
        specFilter <- parseSpeciesFilter(species, wavFiles[w])
        if(!is.null(specFilter)) {
            filtFrom <- specFilter$from
            filtTo <- specFilter$to
        } else {
            filtFrom <- from
            filtTo <- to
        }

        thisWav <- seewave::bwfilter(thisWav, from=filtFrom, to=filtTo)[,1]
        # wavMat[, w] <- thisWav
        wavList[[w]] <- thisWav
    }
    if(length(tooShort) > 0) {
        warning('Plots not created for ', length(tooShort), ' clips that were shorter',
                ' than the desired clip length.', call.=FALSE)
        wavFiles <- wavFiles[-tooShort]
        wavList <- wavList[-tooShort]
    }
    cat('\n')
    cat('Creating plots...\n')
    pb <- txtProgressBar(min=0, max=length(wavFiles), style=3)
    wavEvents <- parseWavClipName(wavFiles, 'Event')
    wavIx <- 1
    # browser()
    for(e in unique(wavEvents)) {
        thisEvent <- wavEvents == e
        thisWavFiles <- wavFiles[thisEvent]
        # thisMat <- wavMat[, thisEvent]
        thisMat <- wavList[thisEvent]
        if(is.null(nrow) &&
           is.null(ncol)) {
            nrow <- ceiling(sqrt(length(thisWavFiles)))
        }
        if(is.null(nrow)) {
            nrow <- ceiling(length(thisWavFiles)/ncol)
        }
        if(is.null(ncol)) {
            ncol <- ceiling(length(thisWavFiles)/nrow)
        }

        # for(i in 1:(ceiling(ncol(thisMat)/(nrow*ncol)))) {
        for(i in 1:(ceiling(length(thisMat)/(nrow*ncol)))) {
            wStart <- 1+ (i-1)*(nrow*ncol)
            # wEnd <- min(ncol(thisMat), nrow*ncol*i)
            wEnd <- min(length(thisMat), nrow*ncol*i)
            # cat('\n',min(ceiling(wEnd/ncol), nrow)*2)
            thisRow <- min(ceiling(wEnd/ncol), nrow)
            png(filename=file.path(outDir, paste0(e,'_WaveformPanels_', wStart, '-', wEnd, '.png')),
                res=300, width=ncol*3, height=thisRow*2, units='in')
            on.exit(dev.off())
            oldMar <- par()$mar
            oldMf <- par()$mfrow
            par(mfrow=c(thisRow, ncol), mar = c(1,1,1,1))
            for(w in wStart:wEnd) {
                setTxtProgressBar(pb, value=wavIx)
                wavIx <- wavIx + 1
                # thisWav <- thisMat[, w]
                thisWav <- thisMat[[w]]
                plot(y=thisWav, x=(1:length(thisWav))/thisSr*1e3, type='l', xlab= 'Time (ms)', yaxt='n')
                if(length(vlines) > 0) {
                    for(v in seq_along(vlines)) {
                        lines(x=rep(vlines[v]*1e3, 2), y=range(thisWav), lty=2)
                    }
                }
                thisDepth <- parseArrayDepth(arrayDepth, thisWavFiles[w], error=FALSE)
                thisDepthErr <- parseArrayDepth(arrayDepth, thisWavFiles[w], error=TRUE)
                if(!is.na(thisDepth) &&
                   length(thisDepth) == 1) {
                    maxEcho <- 2*thisDepth / soundSpeed * 1e3
                    maxErr <- 2*thisDepthErr / soundSpeed * 1e3
                    clickTime <- which.max(abs(thisWav))/thisSr*1e3

                    lines(x=rep(clickTime + maxErr, 2), y=range(thisWav), lty=2, col='dodgerblue')
                    lines(x=rep(clickTime + maxEcho, 2), y=range(thisWav), lty=2, col='blue')

                    # text(x=clickTime + maxEcho, y=max(thisWav, na.rm=TRUE)*.7, col='blue',
                    #      label=paste0('Max Echo (', round(thisDepth, 2), 'm)'), pos=4)
                }

                # maxEcho <- 2*arrayDepth / soundSpeed * 1e3
                # clickTime <- which.max(abs(thisWav))/thisSr*1e3
                # lines(x=rep(clickTime+maxEcho, 2), y=range(thisWav), lty=2, col='blue')
                # text(x=clickTime + maxEcho, y=max(thisWav, na.rm=TRUE)*.7, col='blue', label='Max Echo', pos=4)
                label <- gsub('.*\\.([^\\.]*)\\.wav$', '\\1', thisWavFiles[w])
                evId <- strsplit(parseWavClipName(thisWavFiles[w], 'Event'), '\\.')[[1]][2]
                uid <- parseWavClipName(thisWavFiles[w], 'UID')
                utc <- parseWavClipName(thisWavFiles[w], 'UTC')
                utcMillis <- as.character(round(as.numeric(utc)-floor(as.numeric(utc)), 3))
                utcMillis <- substr(utcMillis, 3, min(5, nchar(utcMillis)))
                # label <- paste0(evId, '.', label)
                evLabel <- paste0(evId, ' UID#', uid)
                timeLabel <- paste0(as.character(utc), '.', utcMillis)
                text(x=300/thisSr*1e3, y = max(thisWav, na.rm=TRUE) * .9, label=evLabel, pos=4)
                text(x=300/thisSr*1e3, y = max(thisWav, na.rm=TRUE) * .75, label=timeLabel, pos=4)
                text(pos=2, x=length(thisWav)/thisSr*1e3, y=min(thisWav, na.rm=TRUE)*.9, label='(ms)', offset=.4)
            }
            par(mfrow=oldMf, mar=oldMar)
            dev.off()
        }
    }
    cat('\n')
    on.exit()
    # wavMat
    wavList
}

# x is wav filename created by createEventClips
# mode is what piece to read from that
parseWavClipName <- function(x, mode=c('UID', 'UTC', 'Event')) {
    if(length(x) > 1) {
        return(sapply(x, function(i) parseWavClipName(i, mode), USE.NAMES = FALSE))
    }
    splitName <- strsplit(basename(x), '\\.')[[1]]
    switch(match.arg(mode),
           UID = {
               gsub('([0-9]*)CH[0-9_]*', '\\1', splitName[length(splitName)-1])
           },
           UTC = {
               utcStr <- splitName[length(splitName)-1]
               milli <- as.numeric(substr(utcStr, stop=nchar(utcStr), start=nchar(utcStr)-2))/1e3
               utcStr <- substr(utcStr, stop=nchar(utcStr)-4, start=nchar(utcStr)-18)
               as.POSIXct(utcStr, format='%Y%m%d_%H%M%S', tz='UTC') + milli
           },
           Event = {
               paste0(gsub('^Event_|^Detection_', '', splitName[1]),
                      '.',
                      splitName[2])
           },
           NA
    )
}

parseArrayDepth <- function(depth, wavName, error=FALSE) {
    if(is.numeric(depth) && length(depth) == 1) {
        return(depth)
    }
    if(is.data.frame(depth)) {
        evName <- parseWavClipName(wavName, 'Event')
        thisUID <- parseWavClipName(wavName, 'UID')
        thisMatch <- depth$UID == thisUID &
            depth$EventId == evName
        thisDepth <- unique(depth$ArrayDepth[thisMatch])
        if(error) {
            return(thisDepth * (1 + unique(depth$ArrayDepthPercentErr[thisMatch]) / 100))
        }
        return(thisDepth)
    }
    if(is.AcousticStudy(depth)) {
        evName <- parseWavClipName(wavName, 'Event')
        thisEv <- depth[evName]
        clickData <- getClickData(thisEv)
        thisUID <- parseWavClipName(wavName, 'UID')
        return(unique(clickData$hpDepth[clickData$UID == thisUID]))
    }
    NA
}

parseSpeciesFilter <- function(x, wavName) {
    # if(is.character(x) && length(x) == 1) {
    #
    # }
    if(is.null(x)) {
        return(NULL)
    }
    if(is.AcousticStudy(x)) {
        evName <- parseWavClipName(wavName, 'Event')
        thisEv <- x[[evName]]
        if(is.null(thisEv)) {
            return(NULL)
        }
        x <- species(thisEv)$id
    }
    # do automatic filter levels by species
    # Cuvier's: 18-50 kHz
    #       True's: 30-70 kHz
    # Gervais'/MmMe: 25-70 kHz
    #       Sowerby's: 50-80 kHz
    # Blainville's: 20-50 kHz
    if(is.na(x)) {
        return(NULL)
    }
    specList <- c('cuviers', 'trues', 'gervais', 'mmme', 'sowerbys', 'blainvilles')
    if(!is.null(x)) {
        x <- match.arg(gsub('[[:punct:]]', '', tolower(x)), choices=specList)
        switch(x,
               cuviers = {
                   # from <- 18e3; to <- 50e3
                   return(list(from = 18e3, to = 50e3))
               },
               trues = {
                   # from <- 30e3; to <- 70e3
                   return(list(from = 30e3, to = 70e3))
               },
               gervais = {
                   # from <- 25e3; to <- 70e3
                   return(list(from = 25e3, to = 70e3))
               },
               mmme = {
                   # from <- 25e3; to <- 70e3
                   return(list(from = 25e3, to = 70e3))
               },
               sowerbys = {
                   # from <- 50e3; to <- 80e3
                   return(list(from = 50e3, to = 80e3))
               },
               blainvilles = {
                   # from <- 20e3; 50e3
                   return(list(from = 20e3, to = 50e3))
               }
               # },
               # warning('No match for species ', species)
        )
    }
    return(NULL)
}

tryMatchNEFSC <- function(x, df, startOnly=TRUE) {
    x <- gsub('[[:punct:]]', '', x)
    x <- tolower(x)
    x <- gsub('^off eff', '', x)
    x <- gsub('^ ', '', x)
    matchStart <- sapply(paste0('^', df$comment), function(c) {
        grepl(c, x)
    })
    if(any(matchStart)) {
        possSpec <- unique(df$id[matchStart])
        if(length(possSpec) == 1) {
            return(possSpec)
        } else {
            return('Multiple')
        }
    }
    if(isFALSE(startOnly)) {
        matchRest <- sapply(df$comment, function(c) {
            grepl(c, x)
        })
        if(!any(matchRest)) {
            return(NA)
        }
        possSpec <- unique(df$id[matchRest])
        if(length(possSpec) == 1) {
            return(possSpec)
        }
        return('Multiple')
    }
    NA
}
# MmMe possible gervais -> gervais currently
# What to do with possible
# Zc a couple times

commentToSpecies <- function(x, validate = FALSE, specMap=NULL) {
    if(is.null(specMap)) {
        specMap <- PAMpal:::ppVars()$specMap
    }
    df <- bind_rows(lapply(events(x), function(e) {
        list(event=id(e),
             comment=ancillary(e)$eventComment)
    }))
    df$species <- unname(sapply(df$comment, function(c) tryMatchNEFSC(c, specMap, startOnly=TRUE)))
    isNA <- is.na(df$species)
    if(any(isNA)) {
        warning(sum(isNA), ' events comments could not be matched to a species,',
                ' these have been set to NA', call.=FALSE)
    }
    if(!all(isNA) &&
       any(df$species[!isNA] == 'Multiple')) {
        warning(sum(df$species[!isNA] == 'Multiple'), ' event comments matched multiple',
                ' species, these have been set to "Multiple"', call.=FALSE)
    }
    if(!isFALSE(validate)) {
        if(is.character(validate)) {
            csvName <- validate
        } else {
            csvName <- 'SpeciesValidation.csv'
        }
        cat('Creating species validation file "', csvName,'"\n', sep='')
        cat('(', sum(isNA) +sum(df$species[!isNA] == 'Multiple'), ' species IDs of NA or "Multiple" to check)\n', sep='')
        write.csv(df, file=csvName, row.names = FALSE)
    }
    df
}
