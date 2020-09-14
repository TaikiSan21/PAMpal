#' @title Create Wav Clips of Events
#'
#' @description Creates an audio clip containing sounds from an event
#'
#' @param event \linkS4class{AcousticEvent} or \linkS4class{AcousticStudy}
#'   object containing events to make wav clips for
#' @param wavFolder the folder where the wav files are located
#' @param buffer amount around each event to also include in the clip, in seconds
#' @param format either \code{'pamguard'} or \code{'soundtrap'} specifying
#'   how the original audio files were created. The names of these files must
#'   not have been changed, the original formatting is used to parse out the
#'   start times of the audio files
#' @param log optional location of SoundTrap XML log files, only relevant for
#'   \code{format='soundtrap'}
#'
#' @return A vector of file names for the wav clips that were successfully
#'   created, any that were not able to be written will be \code{NA}
#'
#' @author Taiki Sakai \email{taiki.sakai@@noaa.gov}
#'
#' @importFrom dplyr bind_rows
#' @importFrom tuneR readWave writeWave MCnames bind
#' @importFrom xml2 read_xml xml_find_all
#'
#' @export
#'
writeEventClips <- function(event, wavFolder=NULL, buffer = 0.1, format=c('pamguard', 'soundtrap'), log=NULL) {
    if(is.null(wavFolder)) {
        # wavFolder <- choose.dir(caption='Select a folder containing your wav files.')
        wavFolder <- tk_choose.dir(caption = 'Select a folder containing your wav files.',
                                   default = getwd())
    }
    if(!dir.exists(wavFolder)) {
        stop('Cannot locate wavFolder.')
    }
    format <- match.arg(format)
    wavs <- list.files(wavFolder, full.names=TRUE)
    if(format == 'soundtrap') {
        if(is.null(log)) {
            # log <- choose.dir(caption='Select a folder of Soundtrap log files (optional)')
            log <- tk_choose.dir(caption = 'Select a folder of SoundTrap log files (optional)',
                                 default = getwd())
        }
        if(dir.exists(log)) {
            stLog <- getSoundtrapLog(log)
        } else {
            stLog <- data.frame(micros=0, sample=1, file=basename(wavs))
        }
    }
    wavMap <- bind_rows(lapply(wavs, function(x) {
        rng <- getWavDate(x, format)
        list(start=rng[1], end=rng[2], file=x, length=as.numeric(difftime(rng[2], rng[1], units='secs')))
    }))
    wavMap <- arrange(wavMap, .data$start)
    wavMap$wavGroup <- 1
    wavMap$timeDiff <- 0
    if(nrow(wavMap) > 1) {
        wg <- 1
        for(i in 2:nrow(wavMap)) {
            wavMap$timeDiff[i] <- as.numeric(difftime(wavMap$start[i], wavMap$end[i-1], units='secs'))
            if(wavMap$timeDiff[i] < 0) {
                wavMap$end[i-1] <- wavMap$start[i]
            }
            if(wavMap$end[i-1] != wavMap$start[i]) {
                wg <- wg + 1
            }
            wavMap$wavGroup[i] <- wg
        }
    }
    # browser()
    if(is.AcousticEvent(event)) {
        event <- list(event)
    }
    if(is.AcousticStudy(event)) {
        event <- events(event)
    }
    allFiles <- vector('character', length = length(event))
    cat('Writing wav files...\n')
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
                warning('Could not find matching wav files for event', event[[i]]@id)
                allFiles[[i]] <- NA
                next
            }
            evRange[1] <- wavMap$start[startIx]
        }
        endIx <- checkIn(evRange[2], wavMap)
        if(is.na(endIx)) {
            endIx <- checkIn(evRange[2] - buffer, wavMap)
            if(is.na(endIx)) {
                warning('Could not find matching wav files for event', event[[i]]@id)
                allFiles[[i]] <- NA
                next
            }
            evRange[2] <- wavMap$end[endIx]
        }
        if(wavMap$wavGroup[startIx] != wavMap$wavGroup[endIx]) {
            warning('Event', event[[i]]@id, 'spanned two non-consecutive wav files, could not create clip.')
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
        writeWave(wavResult, fileName, extensible = FALSE)
        allFiles[[i]] <- fileName
        setTxtProgressBar(pb, value=i)
    }
    isNa <- is.na(allFiles)
    cat('\n', paste0('Wrote ', sum(!isNa), ' wav file(s).'))
    names(allFiles) <- sapply(event, function(x) x@id)
    allFiles
}

checkIn <- function(time, map) {
    possible <- (time >= map$start) & (time <= map$end)
    if(!any(possible)) {
        return(NA)
    }
    which(possible)
}

# wav file name to c(start, end) in posix time
getWavDate <- function(wav, format=c('pamguard', 'soundtrap')) {
    header <- readWave(wav, header = TRUE)
    len <- header$samples / header$sample.rate
    # browser()

    switch(
        format,
        'pamguard' = {
            date <- gsub('.*([0-9]{8}_[0-9]{6}_[0-9]{3})\\.wav$', '\\1', wav)
            posix <- as.POSIXct(substr(date, 1, 15), tz = 'UTC', format = '%Y%m%d_%H%M%S')
            millis <- as.numeric(substr(date, 17, 19)) / 1e3
        },
        'soundtrap' = {
            date <- gsub('.*\\.([0-9]{12})\\.wav$', '\\1', wav)
            posix <- as.POSIXct(date, format = '%y%m%d%H%M%S', tz='UTC')
            millis <- 0
        }
    )
    if(is.na(posix)) {
        stop('Could not convert the name of the wav file to time properly.')
    }
    c(0, len) + posix + millis
}

getSoundtrapLog <- function(x) {
    xFold <- list.files(x, full.names = TRUE, pattern = 'xml')
    missing <- lapply(xFold, function(x) {
        xml <- read_xml(x)
        info <- as.character(xml_find_all(xml, '//@Info'))
        hasSG <- grepl('Sampling Gap', info)
        if(any(hasSG)) {
            return(bind_rows(lapply(info[hasSG], function(i) {
                sg <- as.numeric(strsplit(gsub('.*Sampling Gap ([0-9]*) us at sample ([0-9]*) .*', '\\1_\\2', i), '_')[[1]])
                list(micros=sg[1], sample=sg[2], file = gsub('\\.log\\.xml$', '', basename(x)))
            })))
        }
        data.frame(micros=0, sample=1, file = gsub('\\.log\\.xml$', '', basename(x)), stringsAsFactors = FALSE)
    })
    bind_rows(missing)
}
