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
        wavMap <- dbMap[[d]]
        allFiles <- vector('character', length = length(event))
        cat('Writing wav files for database ', names(dbMap)[d], ' ...\n', sep='')
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
                    next
                }
                evRange[2] <- wavMap$end[endIx]
            }
            if(wavMap$fileGroup[startIx] != wavMap$fileGroup[endIx]) {
                warning('Event ', event[[i]]@id, ' spanned two non-consecutive wav files, could not create clip.')
                allFiles[[i]] <- NA_character_
                setTxtProgressBar(pb, value=i)
                next
            }
            if(any(!file.exists(wavMap$file[startIx:endIx]))) {
                warning('Wav files for event ', event[[i]]@id, ' could not be found on disk.',
                        ' Function "updateFiles" can help relocate files that have moved.')
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
            wavResult <- do.call(bind, wavResult)[, 1:min(2, ncol(wavResult))]
            colnames(wavResult) <- MCnames$name[1:min(2, ncol(wavResult))]
            fileName <- paste0('Event_', event[[i]]@id, '.wav')
            fileName <- paste0(gsub('\\.wav$', '', fileName), '.wav')
            fileName <- file.path(outDir, fileName)
            writeWave(wavResult, fileName, extensible = FALSE)
            allFiles[[i]] <- fileName
            setTxtProgressBar(pb, value=i)
        }
        isNa <- is.na(allFiles)
        cat('\n', paste0('Wrote ', sum(!isNa), ' wav file(s).'))
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
#
# # wav file name to c(start, end) in posix time
# getWavDate <- function(wav, tryFirst=NULL) {
#     header <- readWave(wav, header = TRUE)
#     len <- header$samples / header$sample.rate
#     format <- c(tryFirst, c('pamguard', 'soundtrap', 'sm3'))
#     for(f in format) {
#         switch(
#             f,
#             'pamguard' = {
#                 date <- gsub('.*([0-9]{8}_[0-9]{6}_[0-9]{3})\\.wav$', '\\1', wav)
#                 posix <- as.POSIXct(substr(date, 1, 15), tz = 'UTC', format = '%Y%m%d_%H%M%S')
#                 if(is.na(posix)) next
#                 millis <- as.numeric(substr(date, 17, 19)) / 1e3
#                 if(!is.na(posix)) {
#                     FOUNDFORMAT <<- f
#                     break
#                 }
#             },
#             'soundtrap' = {
#                 date <- gsub('.*\\.([0-9]{12})\\.wav$', '\\1', wav)
#                 posix <- as.POSIXct(date, format = '%y%m%d%H%M%S', tz='UTC')
#                 millis <- 0
#                 if(!is.na(posix)) {
#                     FOUNDFORMAT <<- f
#                     break
#                 }
#             },
#             'sm3' = {
#                 date <- gsub('.*\\_([0-9]{8}_[0-9]{6})\\.wav$', '\\1', wav)
#                 posix <- as.POSIXct(date, format = '%Y%m%d_%H%M%S', tz='UTC')
#                 millis <- 0
#                 if(!is.na(posix)) {
#                     FOUNDFORMAT <<- f
#                     break
#                 }
#             }
#         )
#     }
#
#     if(is.na(posix)) {
#         warning('Could not convert the name of the wav file ',
#                 wav, ' to time properly.', call. = FALSE)
#         return(c(NA, NA))
#     }
#     c(0, len) + posix + millis
# }
#
# getSoundtrapLog <- function(x) {
#     xFold <- list.files(x, full.names = TRUE, pattern = 'xml')
#     missing <- lapply(xFold, function(x) {
#         xml <- read_xml(x)
#         info <- as.character(xml_find_all(xml, '//@Info'))
#         hasSG <- grepl('Sampling Gap', info)
#         if(any(hasSG)) {
#             return(bind_rows(lapply(info[hasSG], function(i) {
#                 sg <- as.numeric(strsplit(gsub('.*Sampling Gap ([0-9]*) us at sample ([0-9]*) .*', '\\1_\\2', i), '_')[[1]])
#                 list(micros=sg[1], sample=sg[2], file = gsub('\\.log\\.xml$', '', basename(x)))
#             })))
#         }
#         data.frame(micros=0, sample=1, file = gsub('\\.log\\.xml$', '', basename(x)), stringsAsFactors = FALSE)
#     })
#     bind_rows(missing)
# }
#
# mapWavFolder <- function(wavFolder=NULL, format=c('pamguard', 'soundtrap', 'sm3'), log=NULL) {
#     if(is.null(wavFolder)) {
#         wavFolder <- tk_choose.dir(caption = 'Select a folder containing your wav files.',
#                                    default = getwd())
#     }
#     if(!dir.exists(wavFolder)) {
#         stop('Cannot locate wavFolder.')
#     }
#     if(length(format) != 1) {
#         fmtChoice <- menu(choices=c('Pamguard', 'SoundTrap', 'SM3'),
#                           title = 'What is the source of your sound files?')
#         if(fmtChoice == 0) {
#             stop('Currently only works with Pamguard, SoundTrap, or SM3 files')
#         }
#         format <- c('pamguard', 'soundtrap', 'sm3')[fmtChoice]
#     }
#     format <- match.arg(format)
#     wavs <- list.files(wavFolder, full.names=TRUE, pattern = '\\.wav$', recursive=TRUE)
#     if(length(wavs) == 0) {
#         stop('No wav files found, please check directory.')
#     }
#     if(format == 'soundtrap') {
#         if(is.null(log)) {
#             # log <- choose.dir(caption='Select a folder of Soundtrap log files (optional)')
#             cat('Select a folder of SoundTrap log files (optional)')
#             log <- tk_choose.dir(caption = 'Select a folder of SoundTrap log files (optional)',
#                                  default = getwd())
#         }
#         if(dir.exists(log)) {
#             stLog <- getSoundtrapLog(log)
#             stLog <- group_by(stLog, file) %>%
#                 summarise(gap = sum(micros)/1e6) %>%
#                 ungroup()
#         } else {
#             stLog <- data.frame(gap=0, sample=1, file=basename(wavs))
#         }
#     }
#     FOUNDFORMAT <- NULL
#     wavMap <- bind_rows(lapply(wavs, function(x) {
#         rng <- getWavDate(x, FOUNDFORMAT)
#         if(any(is.na(rng))) {
#             return(NULL)
#         }
#         if(FOUNDFORMAT == 'soundtrap') {
#             hasLog <- which(gsub('\\.wav$', '', basename(x)) == stLog$file)
#             if(length(hasLog) == 1) {
#                 rng[2] <- rng[2] + stLog$gap[hasLog]
#             }
#         }
#         list(start=rng[1], end=rng[2], file=x, length=as.numeric(difftime(rng[2], rng[1], units='secs')))
#     }))
#     wavMap <- arrange(wavMap, .data$start)
#     wavMap
# }
