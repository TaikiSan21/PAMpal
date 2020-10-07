#' @title Create Wav Clips of Events
#'
#' @description Creates an audio clip containing sounds from an event
#'
#' @param x \linkS4class{AcousticStudy} object containing events to make wav clips for
#' @param buffer amount around each event to also include in the clip, in seconds
#' @param outDir directory to write clips to, defaults to current directory
#'
#' @return A vector of file names for the wav clips that were successfully
#'   created, any that were not able to be written will be \code{NA}
#'
#' @examples
#'
#' data(exStudy)
#' recs <- system.file('extdata', 'Recordings', package='PAMpal')
#' exStudy <- addRecordings(exStudy, folder=recs, log=FALSE, progress=FALSE)
#' tdir <- tempdir()
#' wavs <- writeEventClips(exStudy, outDir=tdir)
#' file.remove(wavs)
#'
#' @author Taiki Sakai \email{taiki.sakai@@noaa.gov}
#'
#' @importFrom dplyr bind_rows arrange group_by summarise ungroup
#' @importFrom tuneR readWave writeWave MCnames bind
#' @importFrom xml2 read_xml xml_find_all
#'
#' @export
#'
writeEventClips <- function(x, buffer = 0.1, outDir='.') {
    if(is.null(files(x)$recordings)) {
        stop('No recording files found, use function "addRecordings" first.')
    }
    if(!dir.exists(outDir)) dir.create(outDir)
    evDbs <- sapply(events(x), function(e) files(e)$db)
    dbMap <- split(files(x)$recordings, files(x)$recordings$db)
    for(d in seq_along(dbMap)) {
        event <- events(x)[which(evDbs == names(dbMap)[d])]
        if(length(event) == 0) next
        wavMap <- dbMap[[d]]
        allFiles <- vector('character', length = length(event))
        cat('Writing wav files for database ', basename(names(dbMap)[d]),
            '(', d, ' of ', length(dbMap),') ...\n', sep='')
        pb <- txtProgressBar(min=0, max=length(event), style = 3)
        for(i in seq_along(event)) {
            # browser()
            evRange <- getEventTime(event[[i]]) + buffer * c(-1, 1)
            # for start and end check if in range. if we buffered, try undoing that first.
            # so like if buffer put us before first file, start and beginning of first file instead.
            startIx <- checkIn(evRange[1], wavMap)
            if(is.na(startIx)) {
                startIx <- checkIn(evRange[1] + buffer, wavMap)
                if(is.na(startIx)) {
                    warning('Could not find matching wav files for event ', event[[i]]@id)
                    allFiles[[i]] <- NA
                    setTxtProgressBar(pb, value=i)
                    next
                }
                evRange[1] <- wavMap$start[startIx]
            }
            endIx <- checkIn(evRange[2], wavMap)
            if(is.na(endIx)) {
                endIx <- checkIn(evRange[2] - buffer, wavMap)
                if(is.na(endIx)) {
                    warning('Could not find matching wav files for event ', event[[i]]@id)
                    allFiles[[i]] <- NA
                    setTxtProgressBar(pb, value=i)
                    next
                }
                evRange[2] <- wavMap$end[endIx]
            }
            if(wavMap$fileGroup[startIx] != wavMap$fileGroup[endIx]) {
                warning('Event ', event[[i]]@id,
                        ' spanned two non-consecutive wav files, could not create clip.', call.=FALSE)
                allFiles[[i]] <- NA_character_
                setTxtProgressBar(pb, value=i)
                next
            }
            if(any(!file.exists(wavMap$file[startIx:endIx]))) {
                warning('Wav files for event ', event[[i]]@id, ' could not be found on disk.',
                        ' Function "updateFiles" can help relocate files that have moved.', call. = FALSE)
                allFiles[[i]] <- NA_character_
                setTxtProgressBar(pb, value=i)
                next
            }
            startTime <- as.numeric(difftime(evRange[1], wavMap$start[startIx], units='secs'))
            endTime <- as.numeric(difftime(evRange[2], wavMap$start[endIx], units='secs'))
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
            fileName <- paste0('Event_', event[[i]]@id, '.wav')
            fileName <- paste0(gsub('\\.wav$', '', fileName), '.wav')
            fileName <- file.path(outDir, fileName)
            writeWave(wavResult, fileName, extensible = FALSE)
            allFiles[[i]] <- fileName
            setTxtProgressBar(pb, value=i)
        }
        isNa <- is.na(allFiles)
        cat('\n', paste0('Wrote ', sum(!isNa), ' wav file(s).\n'))
        names(allFiles) <- sapply(event, function(x) x@id)
        allFiles
    }
    invisible(allFiles)
}

checkIn <- function(time, map) {
    possible <- (time >= map$start) & (time <= map$end)
    if(!any(possible)) {
        return(NA)
    }
    which(possible)
}
