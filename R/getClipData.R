#' @title Get Wav Clips of Data
#'
#' @description Reads audio clips containing sounds from events or detections
#'
#' @param x \linkS4class{AcousticStudy} object containing data to read wav clips for
#' @param buffer amount before and after each event to also include in the clip, in seconds.
#'   Can either be a vector of length two specifying how much to buffer before and after
#'   (first number should be negative), or a single value if the buffer amount should be identical.
#' @param mode either \code{'event'} or \code{'detection'} specifying whether to create
#'   wav clips of entire events or individual detections
#' @param channel channel(s) of clips to write
#' @param useSample logical flag to use startSample information in binaries instead of UTC
#'   time for start of detections. This can be slightly more accurate (~1ms) but will take
#'   longer
#' @param progress logical flag to show progress bar
#' @param verbose logical flag to show summary messages
#' @param FUN optional function to apply to wav clips
#' @param \dots optional arguments to pass to \code{FUN}
#'
#' @return A named list of wav clips
#'
#' @examples
#'
#' data(exStudy)
#' recs <- system.file('extdata', 'Recordings', package='PAMpal')
#' exStudy <- addRecordings(exStudy, folder=recs, log=FALSE, progress=FALSE)
#' \dontrun{
#' # not running so that no wav clips are written to disk
#' wavs <- getClipData(exStudy, mode='event')
#' }
#'
#' @author Taiki Sakai \email{taiki.sakai@@noaa.gov}
#'
#' @importFrom dplyr bind_rows arrange group_by summarise ungroup
#' @importFrom tuneR readWave writeWave MCnames bind
#'
#' @export
#'
getClipData <- function(x, buffer = c(0, 0.1), mode=c('event', 'detection'),
                        channel = 1, useSample=FALSE, progress=TRUE, verbose=TRUE, FUN=NULL, ...) {
    if(!is.AcousticStudy(x)) {
        stop('"x" must be an AcousticStudy object.')
    }
    fileExists <- checkRecordings(x)
    mode <- match.arg(mode)

    evDbs <- sapply(events(x), function(e) files(e)$db)
    dbMap <- split(files(x)$recordings, files(x)$recordings$db)

    if(length(buffer) == 1) {
        buffer <- buffer * c(-1, 1)
    }
    buffer <- abs(buffer) * c(-1, 1)

    if(is.null(FUN)) {
        FUN <- function(wav, ...) wav
    }
    # each database can have a different set of assigned recordings,
    # so we break up by DB

    # setup warning storage
    noMatch <- character(0)
    multMatch <- character(0)
    nonConsec <- character(0)
    fileDNE <- character(0)
    noChan <- character(0)
    result <- character(0)
    on.exit({
        if(length(noMatch) > 0) {
            warning('Could not find matching wav files for ', mode, ' ', printN(noMatch, 6), call.=FALSE)
        }
        if(length(multMatch) > 0) {
            warning(oneUpper(mode), ' ', printN(multMatch, 6),
                    ' matched more than 1 possible wav file, could not get clip.', call.=FALSE)
        }
        if(length(nonConsec) > 0) {
            warning(oneUpper(mode), ' ', printN(nonConsec, 6),
                    ' spanned two non-consecutive wav files, could not get clip.', call.=FALSE)
        }
        if(length(fileDNE) > 0) {
            warning('Wav files for ', mode, ' ', printN(fileDNE, 6), ' could not be found on disk.',
                    ' Function "updateFiles" can help relocate files that have moved.', call. = FALSE)
        }
        if(length(noChan) > 0) {
            warning('Wav files for ', mode, ' ', printN(noChan, 6),
                    ' did not have the desired channels.', call.=FALSE)
        }
    })
    # one DB at a time
    for(d in seq_along(dbMap)) {
        thisDbData <- x[which(evDbs == names(dbMap)[d])]
        if(length(events(thisDbData)) == 0) next
        wavMap <- dbMap[[d]]
        allTimes <- getTimeRange(thisDbData, mode=mode, sample=useSample)
        allResult <- vector('list', length = length(allTimes))
        names(allResult) <- names(allTimes)
        if(progress) {
            cat('Processing wav files for database ', basename(names(dbMap)[d]),
                '(', d, ' of ', length(dbMap),') ...\n', sep='')
            pb <- txtProgressBar(min=0, max=length(allResult), style = 3)
        }
        for(i in seq_along(allResult)) {
            timeRange <- c(allTimes[[i]]$start, allTimes[[i]]$end) + buffer
            # for start and end check if in range. if we buffered, try undoing that first.
            # so like if buffer put us before first file, start and beginning of first file instead.
            startIx <- checkIn(timeRange[1], wavMap)
            if(length(startIx) > 1) {
                multMatch <- c(multMatch, names(allResult)[i])
                if(progress) {
                    setTxtProgressBar(pb, value=i)
                }
                next
            }
            if(is.na(startIx)) {
                startIx <- checkIn(timeRange[1] - buffer[1], wavMap)
                if(is.na(startIx)) {
                    noMatch <- c(noMatch, names(allResult)[i])
                    if(progress) {
                        setTxtProgressBar(pb, value=i)
                    }
                    next
                }
                timeRange[1] <- wavMap$start[startIx]
            }
            endIx <- checkIn(timeRange[2], wavMap)
            if(length(endIx) > 1) {
                multMatch <- c(multMatch, names(allResult)[i])
                if(progress) {
                    setTxtProgressBar(pb, value=i)
                }
                next
            }
            if(is.na(endIx)) {
                endIx <- checkIn(timeRange[2] - buffer[2], wavMap)
                if(is.na(endIx)) {
                    noMatch <- c(noMatch, names(allResult)[i])
                    if(progress) {
                        setTxtProgressBar(pb, value=i)
                    }
                    next
                }
                timeRange[2] <- wavMap$end[endIx]
            }
            if(wavMap$fileGroup[startIx] != wavMap$fileGroup[endIx]) {
                nonConsec <- c(nonConsec, names(allResult)[i])
                if(progress) {
                    setTxtProgressBar(pb, value=i)
                }
                next
            }
            if(any(!file.exists(wavMap$file[startIx:endIx]))) {
                fileDNE <- c(fileDNE, names(allResult)[i])
                if(progress) {
                    setTxtProgressBar(pb, value=i)
                }
                next
            }
            startTime <- as.numeric(difftime(timeRange[1], wavMap$start[startIx], units='secs'))
            endTime <- as.numeric(difftime(timeRange[2], wavMap$start[endIx], units='secs'))
            wavResult <- vector('list', length = endIx)
            for(w in startIx:endIx) {
                readStart <- 0
                readEnd <- Inf
                if(w == startIx) {
                    readStart <- startTime
                }
                if(w == endIx) {
                    readEnd <- endTime
                }
                wavResult[[w]] <- readWave(wavMap$file[w], from = readStart, to = readEnd, units = 'seconds', toWaveMC = TRUE)
            }

            wavResult <- wavResult[!sapply(wavResult, is.null)]
            wavResult <- do.call(bind, wavResult) # [, 1:min(2, ncol(wavResult))]
            chanIn <- channel <= ncol(wavResult@.Data)
            if(!any(chanIn)) {
                noChan <- c(noChan, names(allResult)[i])
                if(progress) {
                    setTxtProgressBar(pb, value = i)
                }
                next
            }
            if(!all(chanIn)) {
                noChan <- c(noChan, names(allResult)[i])
            }
            wavResult <- wavResult[, channel[chanIn]]
            colnames(wavResult) <- MCnames$name[1:ncol(wavResult)]
            allResult[[i]] <- FUN(wavResult, name=names(allResult)[i], time=timeRange, channel=channel[chanIn], mode=mode, ...)

            if(progress) {
                setTxtProgressBar(pb, value=i)
            }
        }
        isNull <- sapply(allResult, is.null)
        if(verbose) {
            cat('\n', paste0('Processed ', sum(!isNull), ' wav file(s).\n'))
        }
        result <- c(result, allResult)
    }
    invisible(result)
}

checkIn <- function(time, map) {
    time <- as.numeric(time)
    possible <- (time >= map$numStart) & (time <= map$numEnd)
    if(!any(possible)) {
        return(NA)
    }
    which(possible)
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
        dets <- distinct(
            bind_rows(lapply(detectors(e), function(d) {
                d[, c('UID', 'UTC'), drop = FALSE]
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
                        return(list(start=thisDate, end=thisDate))
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
                result <- lapply(dets$UTC, function(d) {
                    list(start = d, end = d)
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