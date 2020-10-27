#' @title Create Wav Clips of Data
#'
#' @description Creates audio clips containing sounds from events or detections
#'
#' @param x \linkS4class{AcousticStudy} object containing data to make wav clips for
#' @param buffer amount before and after each event to also include in the clip, in seconds.
#'   Can either be a vector of length two specifying how much to buffer before and after
#'   (first number should be negative), or a single value if the buffer amount should be identical.
#' @param outDir directory to write clips to, defaults to current directory
#' @param mode either \code{'event'} or \code{'detection'} specifying whether to create
#'   wav clips of entire events or individual detections
#' @param progress logical flag to show progress bar
#'
#' @return A vector of file names for the wav clips that were successfully
#'   created, any that were not able to be written will be \code{NA}. Note
#'   that currently this can only write clips with up to 2 channels.
#'
#' @examples
#'
#' data(exStudy)
#' recs <- system.file('extdata', 'Recordings', package='PAMpal')
#' exStudy <- addRecordings(exStudy, folder=recs, log=FALSE, progress=FALSE)
#' \dontrun{
#' # not running so that no wav clips are written to disk
#' wavs <- writeEventClips(exStudy, outDir='WavFolder')
#' }
#'
#' @author Taiki Sakai \email{taiki.sakai@@noaa.gov}
#'
#' @importFrom dplyr bind_rows arrange group_by summarise ungroup
#' @importFrom tuneR readWave writeWave MCnames bind
#' @importFrom xml2 read_xml xml_find_all
#'
#' @export
#'
writeEventClips <- function(x, buffer = c(-0.1, 0.1), outDir='.', mode=c('event', 'detection'), progress=TRUE) {
    if(!is.AcousticStudy(x)) {
        stop('"x" must be an AcousticStudy object.')
    }
    if(is.null(files(x)$recordings)) {
        stop('No recording files found, use function "addRecordings" first.')
    }
    if(!dir.exists(outDir)) dir.create(outDir)
    evDbs <- sapply(events(x), function(e) files(e)$db)
    files(x)$recordings$numStart <- as.numeric(files(x)$recordings$start)
    files(x)$recordings$numEnd <- as.numeric(files(x)$recordings$end)
    dbMap <- split(files(x)$recordings, files(x)$recordings$db)

    if(length(buffer) == 1) {
        buffer <- buffer * c(-1, 1)
    }
    buffer <- abs(buffer) * c(-1, 1)
    # each database can have a different set of assigned recordings,
    # so we break up by DB
    for(d in seq_along(dbMap)) {
        thisDbData <- x[which(evDbs == names(dbMap)[d])]
        if(length(events(thisDbData)) == 0) next
        wavMap <- dbMap[[d]]
        allTimes <- getTimeRange(thisDbData, mode=mode)
        allFiles <- vector('character', length = length(allTimes))
        names(allFiles) <- names(allTimes)
        if(progress) {
            cat('Writing wav files for database ', basename(names(dbMap)[d]),
                '(', d, ' of ', length(dbMap),') ...\n', sep='')
            pb <- txtProgressBar(min=0, max=length(allFiles), style = 3)
        }
        for(i in seq_along(allFiles)) {
            # browser()
            timeRange <- allTimes[[i]] + buffer
            # for start and end check if in range. if we buffered, try undoing that first.
            # so like if buffer put us before first file, start and beginning of first file instead.
            startIx <- checkIn(timeRange[1], wavMap)
            if(is.na(startIx)) {
                startIx <- checkIn(timeRange[1] + buffer[1], wavMap)
                if(is.na(startIx)) {
                    warning('Could not find matching wav files for ', mode, names(allFiles)[i])
                    allFiles[[i]] <- NA
                    if(progress) {
                        setTxtProgressBar(pb, value=i)
                    }
                    next
                }
                timeRange[1] <- wavMap$start[startIx]
            }
            endIx <- checkIn(timeRange[2], wavMap)
            if(is.na(endIx)) {
                endIx <- checkIn(timeRange[2] - buffer[2], wavMap)
                if(is.na(endIx)) {
                    warning('Could not find matching wav files for ', mode, names(allFiles)[i])
                    allFiles[[i]] <- NA
                    if(progress) {
                        setTxtProgressBar(pb, value=i)
                    }
                    next
                }
                timeRange[2] <- wavMap$end[endIx]
            }
            if(wavMap$fileGroup[startIx] != wavMap$fileGroup[endIx]) {
                warning(oneUpper(mode), names(allFiles)[i],
                        ' spanned two non-consecutive wav files, could not create clip.', call.=FALSE)
                allFiles[[i]] <- NA_character_
                if(progress) {
                    setTxtProgressBar(pb, value=i)
                }
                next
            }
            if(any(!file.exists(wavMap$file[startIx:endIx]))) {
                warning('Wav files for ', mode, names(allFiles)[i], ' could not be found on disk.',
                        ' Function "updateFiles" can help relocate files that have moved.', call. = FALSE)
                allFiles[[i]] <- NA_character_
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
            wavResult <- wavResult[, 1:min(2, ncol(wavResult))]
            colnames(wavResult) <- MCnames$name[1:min(2, ncol(wavResult))]
            fileName <- paste0(oneUpper(mode), '_', names(allFiles)[i], '.wav')
            fileName <- paste0(gsub('\\.wav$', '', fileName), '.wav')
            fileName <- file.path(outDir, fileName)
            writeWave(wavResult, fileName, extensible = FALSE)
            allFiles[[i]] <- fileName
            if(progress) {
                setTxtProgressBar(pb, value=i)
            }
        }
        isNa <- is.na(allFiles)
        cat('\n', paste0('Wrote ', sum(!isNa), ' wav file(s).\n'))
        # names(allFiles) <- sapply(event, function(x) x@id)
        allFiles
    }
    invisible(allFiles)
}

checkIn <- function(time, map) {
    time <- as.numeric(time)
    possible <- (time >= map$numStart) & (time <= map$numEnd)
    if(!any(possible)) {
        return(NA)
    }
    which(possible)
}

# named vector for AcEv, or named list of named vectors for AcSt
getTimeRange <- function(x, mode=c('event', 'detection')) {
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
            return(c(start=min(dets$UTC), end=max(dets$UTC)))
        }
        if(mode == 'detection') {
            result <- lapply(dets$UTC, function(d) {
                c(start = d, end = d)
            })
            names(result) <- dets$UID
            return(result)
        }
    })
    if(mode == 'detection') {
        allDets <- unlist(allDets, recursive=FALSE)
    }
    allDets
}

oneUpper <- function(x) {
    paste0(toupper(substr(x, 1, 1)),
           tolower(substr(x, 2, nchar(x))))
}
