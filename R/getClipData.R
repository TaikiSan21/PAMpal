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
#' @param fixLength logical flag to fix the output clip length to a constant value. If
#'   \code{TRUE}, then output clip length is entirely determined by the buffer value, as
#'   if the detection or event had zero length. E.g. \code{buffer=c(-2,1)} will produce clips
#'   3 seconds long, starting 2 seconds before the detection/event start time.
#' @param fillZeroes logical flag to fill gaps in non-consecutive clips with zeroes. If
#'   \code{FALSE}, will give warnings when attempting to retrieve clips spanning non-consecutive
#'   files and return no clip data
#' @param progress logical flag to show progress bar
#' @param verbose logical flag to show summary messages
#' @param FUN optional function to apply to wav clips. This function takes default inputs \code{wav},
#'   a Wave class object, \code{name} the name of the detection or event, \code{time} the start and end
#'   time of the clip, \code{channel} as above, \code{mode} as above, and additional args \dots
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
#' @importFrom tuneR readWave writeWave MCnames bind nchannel
#'
#' @export
#'
getClipData <- function(x, buffer = c(0, 0.1), mode=c('event', 'detection'),
                        channel = 1, useSample=FALSE, fixLength=FALSE,
                        fillZeroes=TRUE,
                        progress=TRUE, verbose=TRUE, FUN=NULL, ...) {
    if(!is.AcousticStudy(x)) {
        stop('"x" must be an AcousticStudy object.')
    }
    recs <- checkRecordings(x)
    mode <- match.arg(mode)

    evDbs <- sapply(events(x), function(e) basename(files(e)$db))
    dbMap <- split(recs, recs$db)
    names(dbMap) <- basename(names(dbMap))

    if(length(buffer) == 1) {
        buffer <- buffer * c(-1, 1)
    }
    # buffer <- abs(buffer) * c(-1, 1)

    if(is.null(FUN)) {
        FUN <- function(wav, ...) wav
    }
    if(fixLength &&
       buffer[2] - buffer[1] <= 0) {
        stop('For fixed length output buffer[2] - buffer[1] must be positive')
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
            if(fillZeroes) {
                # warning(oneUpper(mode), ' ', printN(nonConsec, 6),
                        # ' spanned non-consecutive wav files, clip has been zero-filled.', call.=FALSE)
            } else {
                warning(oneUpper(mode), ' ', printN(nonConsec, 6),
                        ' spanned two non-consecutive wav files, could not get clip.', call.=FALSE)
            }
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
        thisDbMatch <- which(evDbs == names(dbMap)[d])
        if(length(thisDbMatch) == 0) next
        thisDbData <- x[thisDbMatch]
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
            if(fixLength) {
                allTimes[[i]]$end <- allTimes[[i]]$start
            }
            timeRange <- c(allTimes[[i]]$start, allTimes[[i]]$end) + buffer
            # these are for times where desired start/end is outside of wav
            if(fillZeroes) {
                zeroBuff <- c(0,0)
            }
            # for start and end check if in range. if we buffered, try undoing that first.
            # so like if buffer put us before first file, start and beginning of first file instead.
            startIx <- checkIn(timeRange[1], wavMap)
            if(length(startIx) > 1) {
                if(length(startIx) == 2 &&
                   diff(startIx) == 1 &&
                   wavMap$timeDiff[startIx[2]] < 0) {
                    startIx <- startIx[1]
                } else {
                    multMatch <- c(multMatch, names(allResult)[i])
                    if(progress) {
                        setTxtProgressBar(pb, value=i)
                    }
                    next
                }
            }
            # if cant match, try undoing buffer then match
            if(is.na(startIx)) {
                # startIx <- checkIn(timeRange[1] - buffer[1], wavMap)
                # if(is.na(startIx)) {
                # if we arent filling zeroes we have to actually match
                if(!fillZeroes) {
                    noMatch <- c(noMatch, names(allResult)[i])
                    if(progress) {
                        setTxtProgressBar(pb, value=i)
                    }
                    next
                }
                startIx <- min(which(wavMap$start > timeRange[1]))
                # if fill zeroes, need to pre-pend by amount of zeroes here
                if(fillZeroes) {
                    zeroBuff[1] <- as.numeric(difftime(wavMap$start[startIx], timeRange[1], units='secs'))
                }
                timeRange[1] <- wavMap$start[startIx]
            }
            endIx <- checkIn(timeRange[2], wavMap)
            if(length(endIx) > 1) {
                if(length(endIx) == 2 &&
                   diff(endIx) == 1 &&
                   wavMap$timeDiff[endIx[2]] < 0) {
                    endIx <- endIx[1]
                } else {
                    multMatch <- c(multMatch, names(allResult)[i])
                    if(progress) {
                        setTxtProgressBar(pb, value=i)
                    }
                    next
                }
            }
            # if cant match, try undoing buffer then match
            if(is.na(endIx)) {
                # endIx <- checkIn(timeRange[2] - buffer[2], wavMap)
                # if(is.na(endIx)) {
                if(!fillZeroes) {
                    noMatch <- c(noMatch, names(allResult)[i])
                    if(progress) {
                        setTxtProgressBar(pb, value=i)
                    }
                    next
                }
                # if fill zeroes, backfill by amount of zeroes here
                endIx <- max(which(wavMap$end < timeRange[2]))
                if(fillZeroes) {
                    zeroBuff[2] <- as.numeric(difftime(timeRange[2], wavMap$end[endIx], units='secs'))
                }
                timeRange[2] <- wavMap$end[endIx]
            }
            if(wavMap$fileGroup[startIx] != wavMap$fileGroup[endIx]) {
                nonConsec <- c(nonConsec, names(allResult)[i])
                if(isFALSE(fillZeroes)) {
                    if(progress) {
                        setTxtProgressBar(pb, value=i)
                    }
                    next
                }
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
                wavResult[[w]] <- readWave(wavMap$file[w], from = readStart,
                                           to = readEnd, units = 'seconds', toWaveMC = TRUE)
                thisSr <- wavResult[[w]]@samp.rate
                if(fillZeroes &&
                   w == startIx &&
                   zeroBuff[1] != 0) {
                    thisZeroes <- WaveMC(data=matrix(0,
                                                     ncol=nchannel(wavResult[[w]]),
                                                     nrow=thisSr*zeroBuff[1]),
                                         samp.rate=thisSr, bit=wavResult[[w]]@bit)
                    wavResult[[w]] <- bind(thisZeroes, wavResult[[w]])
                }
                if(fillZeroes &&
                   w == endIx &&
                   zeroBuff[2] != 0) {
                    thisZeroes <- WaveMC(data=matrix(0,
                                                     ncol=nchannel(wavResult[[w]]),
                                                     nrow=thisSr*zeroBuff[2]),
                                         samp.rate=thisSr, bit=wavResult[[w]]@bit)
                    wavResult[[w]] <- bind(wavResult[[w]], thisZeroes)
                }
                if(fillZeroes &&
                   w != endIx &&
                   wavMap$timeDiff[w+1] > 0) {

                    thisZeroes <- WaveMC(data=matrix(0, ncol=nchannel(wavResult[[w]]), nrow=thisSr*wavMap$timeDiff[w+1]),
                                         samp.rate=thisSr, bit=wavResult[[w]]@bit)
                    wavResult[[w]] <- bind(wavResult[[w]], thisZeroes)
                }
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
