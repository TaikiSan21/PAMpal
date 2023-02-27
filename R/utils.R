# random utils

# for converting from database UTC columns that are characters
pgDateToPosix <- function(x) {
    as.POSIXct(as.character(x), format='%Y-%m-%d %H:%M:%OS', tz='UTC')
}

# drop columns with names cols
dropCols <- function(x, cols) {
    ct <- attr(x, 'calltype')
    keepCols <- !(colnames(x) %in% cols)
    x <- x[, keepCols, drop=FALSE]
    attr(x, 'calltype') <- ct
    x
}

# what event is a UID in returns named index
whereUID <- function(study, UID, quiet=FALSE) {
    UID <- as.character(UID)
    where <- sapply(UID, function(u) { #for each uid
        ix <- which(sapply(events(study), function(e) { #go through each event
            any(sapply(detectors(e), function(d) { #and every detector in that event
                u %in% d$UID
            }))
        }))
        if(length(ix) == 0) {
            return(NA)
        }
        ix
    }, USE.NAMES=TRUE, simplify=FALSE)
    whereNA <- is.na(where)
    if(!quiet && any(whereNA)) {
        pamWarning('UID(s) ', paste0(UID[whereNA], collapse=', '),
                   ' not found in any events.')
    }
    where
}

# match SR function
# data needs UTC, thats it
# db is sound acq data, either DF or db
# safe to fail with NULL instead of error
#' @importFrom data.table data.table setkeyv
#'
matchSR <- function(data, db, extraCols = NULL, safe=FALSE, fixNA=TRUE) {
    if(is.character(db)) {
        if(!file.exists(db)) {
            if(safe) return(NULL)
            stop('Database ', db, ' does not exist.')
        }
        con <-dbConnect(db, drv=SQLite())
        on.exit(dbDisconnect(con))
        if(!('Sound_Acquisition' %in% dbListTables(con))) {
            soundAcquisition <- NULL
        } else {
            soundAcquisition <- dbReadTable(con, 'Sound_Acquisition')
            soundAcquisition$UTC <- pgDateToPosix(soundAcquisition$UTC)
        }
    }
    if(is.data.frame(db)) {
        soundAcquisition <- db
    }
    if(inherits(data, 'POSIXct')) {
        data <- data.frame(UTC = data)
    }
    if(!('UTC' %in% colnames(data)) ||
       !inherits(data$UTC, 'POSIXct')) {
        if(safe) return(NULL)
        stop('Data must have a column "UTC" in POSIXct format.')
    }
    if(!is.null(soundAcquisition)) {
        soundAcquisition <- soundAcquisition %>%
            mutate(Status = str_trim(.data$Status),
                   SystemType = str_trim(.data$SystemType)) %>%
            filter(.data$Status=='Start') %>%
            arrange(.data$UTC) %>%
            select(all_of(c('UTC', 'sampleRate', extraCols))) %>%
            distinct() %>%
            data.table()

        # setkeyv(soundAcquisition, 'UTC')

        data <- data.table(data)
        # setkeyv(data, 'UTC')

        # This rolling join rolls to the first time before. Since we filtered to only starts, it goes back
        # to whatever the last Start was.
        data <- soundAcquisition[data, roll = TRUE, on='UTC'] %>%
            data.frame()
        srNa <- which(is.na(data$sampleRate))
    } else {
        data[extraCols] <- NA
        srNa <- rep(TRUE, nrow(data))
    }
    if(fixNA) {
        if(length(srNa) == nrow(data)) {
            cat('\nNo Sample Rate found in SoundAcquisition table. Enter Sample Rate for this data:\n')
            srReplace <- as.integer(readline())
            data$sampleRate[srNa] <- srReplace
        } else if(length(srNa) > 0) {
            # get mode
            mode <- which.max(tabulate(data$sampleRate[-srNa]))
            srChoice <- menu(title=paste0('Could not get Sample Rate for all detections from the "SoundAcquistion" table.',
                                          ' Should missing values be replaced with ', mode, '(value found in table).'),
                             choices = c('Yes', 'No (I want to enter my own SR)'))
            srReplace <- switch(srChoice,
                                '1' = mode,
                                '2' = {
                                    cat('\nWhat Sample Rate should be used?\n')
                                    readline()
                                }, {
                                    if(safe) return(NULL)
                                    stop('Sample Rate required for calculations.')
                                }
            )
            data$sampleRate[srNa] <- srReplace
        }
    }
    data
}

# check if in start/stop interval
# bounds is a single start/stop, sa is sound acq table from db
#' @importFrom tidyr spread
#'
inInterval <- function(bounds, sa) {
    sa <- sa[sa$Status %in% c('Start', 'Stop'), c('UTC', 'Status', 'sampleRate')]
    if(nrow(sa) < 2 ||
       !all(c('Start', 'Stop') %in% sa$Status)) {
        return(FALSE)
    }
    first <- min(which(sa$Status == 'Start'))
    last <- max(which(sa$Status == 'Stop'))
    if(first > last) {
        return(FALSE)
    }
    sa <- sa[first:last,]
    alt <- sa$Status[1:(nrow(sa)-1)] != sa$Status[2:nrow(sa)]
    sa <- sa[c(TRUE, alt), ]
    sa$id <- rep(1:(nrow(sa)/2), each=2)
    sa <- tidyr::spread(sa, 'Status', 'UTC')
    startIn <- (any((bounds[1] >= sa[['Start']]) & (bounds[1] <= sa[['Stop']])))
    endIn <- (any((bounds[2] >= sa[['Start']]) & (bounds[2] <= sa[['Stop']])))
    contain <- (any((bounds[1] <= sa[['Start']]) & (bounds[2] >= sa[['Stop']])))
    startIn || endIn || contain
}

# add list without replacing old one, only replace matching names
safeListAdd <- function(x, value) {
    if(is.null(value)) {
        return(x)
    }
    if(is.list(value) &&
       length(value) == 0) {
        return(x)
    }
    if(!is.list(value) ||
       is.null(names(value))) {
        stop('Can only add named lists ')
    }
    hasName <- names(value) %in% names(x)
    if(any(hasName)) {
        for(n in names(value)[hasName]) {
            x[[n]] <- value[[n]]
        }
    }
    if(any(!hasName)) {
        x <- c(x, value[!hasName])
    }
    x
}

#' @importFrom lubridate int_standardize
#' 
withinLHS <- function(a, int) {
    int <- int_standardize(int)
    as.numeric(a) - as.numeric(int@start) < int@.Data & as.numeric(a) - as.numeric(int@start) >= 0
}

printN <- function(x, n=6, collapse=', ') {
    nItems <- length(x)
    if(nItems == 0) {
        return('')
    }
    if(nItems > n) {
        x <- c(x[1:n], paste0('... (', nItems-n, ' more not shown)'))
    }
    paste0(paste(x, collapse=collapse))
}

getPamFft <- function(data) {
    if(inherits(data, 'PamBinary')) {
        # data$data <- contourToFreq(data$data)
        return(getPamFft(data$data))
    }
    if(length(data) == 0) {
        return(NULL)
    }
    if(!(all(c('sliceData', 'nSlices', 'sampleDuration', 'startSample', 'maxFreq') %in%
             names(data[[1]])))) {
        # stop('Appears data is not a Whistle and Moan Detector binary file.')
        return(NULL)
    }
    tempData <- data[[1]]
    if(tempData$sliceData[[1]]$sliceNumber == 0) {
        tempData <- data[[2]]
    }
    fftHop <- (tempData$startSample + 1)/tempData$sliceData[[1]]$sliceNumber
    fftLen <- tempData$sampleDuration -
        (tempData$sliceData[[tempData$nSlices]]$sliceNumber - tempData$sliceData[[1]]$sliceNumber) * fftHop
    sr <- fftLen * tempData$maxFreq /
        max(unlist(lapply(tempData$sliceData, function(x) x$peakData)))
    list(sr=sr, hop=fftHop, wl=fftLen)
}

ppVars <- function() {
    list(nonModelVars = c('UID', 'Id', 'parentUID', 'sampleRate', 'Channel',
                          'angle', 'angleError', 'peakTime', 'depth', 'sr'),
         # tarMoCols = c(
         #     "TMModelName1", "TMLatitude1", "TMLongitude1", "BeamLatitude1",
         #     "BeamLongitude1", "BeamTime1", "TMSide1", "TMChi21", "TMAIC1", "TMProbability1",
         #     "TMDegsFreedom1", "TMPerpendicularDistance1", "TMPerpendicularDistanceError1", "TMDepth1",
         #     "TMDepthError1","TMHydrophones1","TMComment1","TMError1","TMLatitude2","TMLongitude2",
             # "BeamLatitude2","BeamLongitude2","BeamTime2","TMSide2", "TMChi22","TMAIC2",
             # "TMProbability2", "TMDegsFreedom2", "TMPerpendicularDistance2", "TMPerpendicularDistanceError2",
             # "TMDepth2" ,"TMDepthError2", "TMHydrophones2","TMError2","TMComment2"),
         tarMoCols = c('TMModelName1', 'TMLatitude1', 'TMLongitude1', 'TMPerpendicularDistance1',
                        'TMPerpendicularDistanceError1', 'TMDepth1', 'TMDepthError1'),
         locCols = c('locName', 'locLat', 'locLong', 'perpDist', 'perpDistErr', 'locDepth', 'depthErr'),
         bftHeight = data.table(bftMax=c(1, 2.4, 2.9, 3.4, 3.9, 4.4, 4.9, 5.4, 6, 12),
                                waveHeight=c(0, 0.05, 0.1, 0.2, 0.5, 0.6, 1.25, 1.3, 2.5, 2.5),
                                key='bftMax'),
         specMap = data.frame(comment = c('cuviers', 'cuvier', 'gervais', 'gervai',
                                          'sowerbys', 'sowerby', 'trues', 'true',
                                          'blainvilles', 'blainville',
                                          'unid mesoplodon', 'mmme', 'undi mesoplodon'),
                              id = c('Cuviers', 'Cuviers', 'Gervais', 'Gervais',
                                     'Sowerbys', 'Sowerbys', 'Trues', 'Trues',
                                     'Blainvilles', 'Blainvilles',
                                     'MmMe', 'MmMe', 'MmMe')),
         dglCols = c('Id', 'UID', 'UTC', 'UTCMilliseconds', 'PCLocalTime', 'PCTime',
                     'ChannelBitmap', 'SequenceBitmap', 'EndTime', 'DataCount'),
         binPattern = '(Clicks|WhistlesMoans|GPL).*pgdf$'
    )
}

# logic to see whether Type or Name is holding wav files
findWavCol <- function(sa) {
    st <- sa$SystemType
    sn <- sa$SystemName
    wavCol <- NA
    if((all(is.na(sn)) || is.null(sn)) &&
       !all(grepl('^Audio ', st))) {
        wavCol <- 'SystemType'
    }
    if(!(all(is.na(sn)) || is.null(sn)) &&
       !(all(grepl('Audio ', sn)))) {
        wavCol <- 'SystemName'
    }
    wavCol
}


getSr <- function(x, type=c('click', 'whistle', 'cepstrum'), name=NULL, data=NULL) {
    # need to vectorize better this will do matchSR over and over again FML
    # if(length(name) > 1) {
    #     return(sapply(name, function(n) {
    #         getSr(x, type, n, UTC)
    #     }))
    # }
    if(!is.AcousticStudy(x) &&
       !is.AcousticEvent(x)) {
        return(NULL)
    }
    # type <- match.arg(type)
    if(length(name) != length(data$UTC)) {
        if(length(name) == 1 &&
           length(data$UTC) > 1) {
            name <- rep(name, length(data$UTC))
        }
        # if(length(data$UTC) == 1 &&
        #    length(name) > 1) {
        #     UTC <- rep(UTC, length(name))
        # }
    }
    srOut <- rep(NA, max(1, length(name)))
    if(is.AcousticEvent(x)) {
        if(length(settings(x)$sr) == 1) {
            return(rep(settings(x)$sr, length(srOut)))
        }
        if(!is.null(data)) {
            # srDf <- bind_rows(lapply(files(x)$db, function(d) {
            #     utcDf <- data.frame(ix=1:length(UTC), UTC=UTC)
            #     utcDf <- matchSR(utcDf, d, safe=TRUE, fixNA=FALSE)
            #     arrange(utcDf, ix)[c('UTC', 'sampleRate')]
            # }))
            data$ix <- 1:nrow(data)
            srDf <- bind_rows(lapply(split(data, data$db), function(d) {
                utcDf <- matchSR(d, d$db[1], safe=TRUE, fixNA=FALSE)
                utcDf
            }))
            srDf <- arrange(srDf, .data$ix)
            srOut <- srDf$sampleRate
            if(!all(is.na(srOut))) {
                return(srOut)
            }
        }
        return(NULL)
    }
    for(i in seq_along(srOut)) {
        srOut[i] <- doOneSr(x, type, name[i])
    }

    srNa <- is.na(srOut)
    if(!any(srNa)) {
        return(srOut)
    }
    # MAKE DOONESR WORK HERE THEN LEFTOVERS GET FIXED WITH UTC FORM DB
    # PORBABLY NEED TO RETURN NA INSTEAD OF NULL, NEED TO CHECK OTHER
    # FUNS THAT USE GETSR TO SEE WHAT THEY CHECK ON FAILURE
    # if we are trying to match by times, do it from the database
    if(!is.null(data)) {
        # srDf <- bind_rows(lapply(files(x)$db, function(d) {
        #     utcDf <- data.frame(ix=1:length(UTC[srNa]), UTC=UTC[srNa])
        #     utcDf <- matchSR(utcDf, d, safe=TRUE, fixNA=FALSE)
        #     arrange(utcDf, ix)[c('UTC', 'sampleRate')]
        #     # matchSR(UTC[srNa], d, safe=TRUE, fixNA=FALSE)
        # }))
        data$ix <- 1:nrow(data)
        srDf <- bind_rows(lapply(split(data, data$db), function(d) {
            utcDf <- matchSR(d[srNa, ], d$db[srNa][1], safe=TRUE, fixNA=FALSE)
            utcDf
        }))
        srDf <- arrange(srDf, .data$ix)
        srOut[srNa] <- srDf$sampleRate
    }
    # giv eup
    if(!all(is.na(srOut))) {
        return(srOut)
    }
    NULL
}

doOneSr <- function(x, type=c('click', 'whistle', 'cepstrum'), name=NULL) {
    detSets <- settings(x)$detectors
    # fixing click detector names to match settings, and filter down if matching name
    if(!is.null(name) &&
       !is.null(detSets)) {
        name <- gsub('_[0-9]{0,3}$', '', name)
        detSets <- detSets[names(detSets) == name]
    }
    # if we have settings, see if any have matching type and one answer
    if(!is.null(detSets) && length(detSets) > 0) {
        whichThisType <- sapply(detSets, function(d) {
            d$type %in% type
        })
        if(any(whichThisType)) {
            possSr <- unique(sapply(detSets[whichThisType], function(d) {
                d$sr
            }))
            if(length(possSr) == 1) {
                return(possSr)
            }
        }
    }
    # try to grab sr_hz param if its clicks
    if('click' %in% type) {
        clickFuns <- pps(x)@functions$ClickDetector
        srSettings <- unique(unlist(lapply(clickFuns, function(c) {
            formals(c)[['sr_hz']]
        })))
        if(length(srSettings) == 1 &&
           is.numeric(srSettings)) {
            return(srSettings)
        }
    }
    # see if theres a single SR from audio settings
    possSr <- unique(sapply(events(x), function(e) {
        settings(e)$sr
    }))
    if(length(possSr) == 1) {
        return(possSr)
    }
    NA
}

# returns named vector for AcEv, or named list of named vectors for AcSt
getTimeRange <- function(x, mode=c('event', 'detection'), sample=FALSE) {
    mode <- match.arg(mode)
    # if(is.AcousticStudy(x)) {
    #     return(lapply(events(x), function(e) {
    #         getTimeRange(e, mode)
    #     }))
    # }
    allDets <- lapply(events(x), function(e) {
        if(length(detectors(e)) == 0) {
            return(NULL)
        }
        dets <- distinct(
            bind_rows(lapply(detectors(e), function(d) {
                if(is.null(d) ||
                   nrow(d) == 0) {
                    return(NULL)
                }
                out <- d[, c('UID', 'UTC', 'duration'), drop = FALSE]
                if('duration' %in% colnames(d)) {
                    switch(attr(d, 'calltype'),
                           'whistle' = out$duration <- d$duration,
                           # 'click' = out$duration <- d$duration / 1e6,
                           'click' = out$duration <- 0, # duration is not reliable for clicks
                           'cepstrum' = out$duration <- d$duration
                    )
                } else {
                    out$duration <- 0
                }
                out
            }))
        )
        if(mode == 'event') {
            if(sample) {
                minUID <- dets$UID[which.min(dets$UTC)[1]]
                maxUID <- dets$UID[which.max(dets$UTC)[1]]
                minUTC <- min(dets$UTC)
                maxUTC <- max(dets$UTC)
                recMap <- files(x)$recordings
                minIx <- checkIn(minUTC, recMap)
                if(is.na(minIx) ||
                   (length(minIx) != 1) ||
                   is.na(recMap$startSample[minIx])) {
                    evResult <- list(start = minUTC)
                } else {
                    binMin <- getBinaryData(x, minUID)[[1]]
                    binSr <- ifelse(is.na(binMin$sr), recMap$sr[minIx], binMin$sr)
                    evResult <- list(start = recMap$start[minIx] +
                                         binMin$startSample/binSr - recMap$startSample[minIx]/recMap$sr[minIx])
                }
                maxIx <- checkIn(maxUTC, recMap)
                if(is.na(maxIx) ||
                   (length(maxIx) != 1) ||
                   is.na(recMap$startSample[maxIx])) {
                    evResult$end <- maxUTC
                } else {
                    binMax <- getBinaryData(x, maxUID)[[1]]
                    binSr <- ifelse(is.na(binMax$sr), recMap$sr[maxIx], binMax$sr)
                    evResult$end <- recMap$start[maxIx] +
                        binMax$startSample / binSr - recMap$startSample[maxIx] / recMap$sr[maxIx]
                }
                return(evResult)
            }
            return(list(start=min(dets$UTC), end=max(dets$UTC)))
        }
        if(mode == 'detection') {
            if(sample) {
                recMap <- files(x)$recordings
                result <- lapply(getBinaryData(x, dets$UID), function(b) {
                    thisDate <- b$date
                    wavIx <- checkIn(thisDate, recMap)
                    if(is.na(wavIx) ||
                       (length(wavIx) != 1) ||
                       is.na(recMap$startSample[wavIx])) {
                        return(list(start=thisDate, end=thisDate + dets$duration[dets$UID == b$UID][1]))
                    }
                    binSr <- ifelse(is.na(b$sr), recMap$sr[wavIx], b$sr)
                    out <- list(start=recMap$start[wavIx] +
                                    b$startSample / binSr - recMap$startSample[wavIx] / recMap$sr[wavIx])
                    out$end <- out$start
                    if('sampleDuration' %in% names(b)) {
                        out$end <- out$end + b$sampleDuration / recMap$sr[wavIx]
                    }
                    out
                })

            } else {
                # result <- lapply(dets$UTC, function(d) {
                #     list(start = d, end = d)
                # })
                result <- lapply(1:nrow(dets), function(d) {
                    list(start = dets$UTC[d], end = dets$UTC[d] + dets$duration[d])
                })
                names(result) <- dets$UID
            }
            return(result)
        }
    })
    if(mode == 'detection') {
        allDets <- unlist(allDets, recursive=FALSE)
    }
    allDets
}

checkSameDetections <- function(x, y) {
    xDet <- getDetectorData(x)
    yDet <- getDetectorData(y)
    if(!all(names(xDet) %in% names(yDet)) ||
       !all(names(yDet) %in% names(xDet))) {
        warning('Different detectors')
        return(FALSE)
    }
    for(d in names(xDet)) {
        if(nrow(xDet[[d]]) != nrow(yDet[[d]])) {
            warning('Different number of detections for detector ', d)
            return(FALSE)
        }
        xy <- nrow(setdiff(xDet[[d]], yDet[[d]]))
        if(xy != 0) {
            warning('xy setdiff is ', xy)
            return(FALSE)
        }
        yx <- nrow(setdiff(yDet[[d]], xDet[[d]]))
        if(yx != 0) {
            warning('yx setdiff is ', yx)
            return(FALSE)
        }
    }
    TRUE
}