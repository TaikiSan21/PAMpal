#' @title Add Recordings to an AcousticStudy Object
#'
#' @description Adds recording files to an AcousticStudy object, runs
#'   interactively to allow users to select files if they are not
#'   provided. No actual
#'   recordings are stored, a dataframe containing information on the
#'   start and end times of the recording files is added to the object.
#'
#' @details \code{mapWavFolder} returns a dataframe of start and end
#'   times
#'
#' @param x a \linkS4class{AcousticStudy} object to add recordings to
#' @param folder a folder of recordings to add. If \code{NULL}, user will be
#'   prompted to select a folder of recordings for each database present in
#'   \code{x}. If a single folder, this will be applied to all databases. If
#'   multiple folders, length must be equal to the number of databases and they
#'   will be applied to each database in the provided order.
#' @param log (optional) log files for SoundTrap recordings. These are used to
#'   adjust apparent lengths of recordings for missing data. If \code{NULL}, user
#'   will be prompted to provide a folder (selecting no folder is a valid option here).
#'   If \code{FALSE} this step will be skipped. If a single folder or multiple folders
#'   will be applied similar to \code{folder}
#' @param progress logical flag to show progress bars
#'
#' @return the same object as \code{x} with recording information added
#'   to the \code{files} slots. The information added is a dataframe containing
#'   the start and end times of recording
#'
#' @author Taiki Sakai \email{taiki.sakai@@noaa.gov}
#'
#' @examples
#'
#' data(exStudy)
#' recs <- system.file('extdata', 'Recordings', package='PAMpal')
#' exStudy <- addRecordings(exStudy, folder=recs, log=FALSE, progress=FALSE)
#' files(exStudy)$recordings
#'
#' @importFrom tcltk tk_choose.dir
#' @export
#'
addRecordings <- function(x, folder=NULL, log=FALSE, progress=TRUE) {
    dbMap <- vector('list', length = length(files(x)$db))
    names(dbMap) <- files(x)$db
    if(is.null(log)) {
        logMsg <- paste0('Are there SoundTrap log files for any of your recordings?')
        logChoice <- menu(choices=c('Yes', 'No'), title = logMsg)
        log <- logChoice == 1
    }
    if(length(log) == 1) {
        log <- rep(log, length(dbMap))
    }
    if(length(log) != length(dbMap)) {
        stop('Number of log folders must either be 1 or equal to number of databases.')
    }
    if(is.null(folder)) {
        folder <- vector('character', length = length(dbMap))
        for(d in seq_along(files(x)$db)) {
            dbMsg <- paste0('Select a recording folder for database ',
                            basename(files(x)$db[d]))
            cat(dbMsg, '\n')
            folder[d] <- tk_choose.dir(caption = dbMsg, default = getwd())
        }
    }
    for(l in seq_along(log)) {
        #optional add log, set to false if no dir selected since not
        # all have to be ST recordings
        # FML this is so janky
        if(isTRUE(log[l]) ||
           log[l] == 'TRUE') {
            logMsg <- paste0('Select a log folder for database ',
                             basename(files(x)$db[l]))
            cat(logMsg, '\n')
            log[l] <- tk_choose.dir(caption = logMsg, default=getwd())
            if(is.na(log[l])) {
                log[l] <- FALSE
            }
        }
    }
    # these should only happen if manually provided, check proper input length
    if(length(folder) == 1) {
        folder <- rep(folder, length(dbMap))
    }

    if(length(folder) != length(dbMap)) {
        stop('Number of folders must either be 1 or equal to the number of databases.')
    }

    logList <- vector('list', length = length(unique(log)))
    names(logList) <- as.character(unique(log))
    for(l in seq_along(logList)) {
        logList[[l]] <- getSoundtrapLog(unique(log)[l])
    }
    isMapped <- NULL
    for(d in seq_along(dbMap)) {
        if(folder[d] %in% isMapped) {
            wavMap <- dbMap[[min(which(folder == folder[d]))]]
        } else {
            wavMap <- mapWavFolder(folder[d], log=logList[[as.character(log[d])]], progress)
            # Try to do start sample
            if(!is.null(wavMap)) {
                if(!file.exists(names(dbMap)[d])) {
                    pamWarning('Database ', names(dbMap)[d], ' could not be found, "startSample"',
                               ' cannot be set.')
                    wavMap$startSample <- 1
                } else {
                    sa <- readSa(names(dbMap)[d])
                    wavCol <- findWavCol(sa)
                    if(!is.na(wavCol)) {
                        wavMap$startSample <- NA
                        for(w in seq_along(wavMap$file)) {
                            thisWav <- grep(substr(basename(wavMap$file[w]), 1, 49), sa[[wavCol]])
                            if(length(thisWav) == 0) next
                            thisSa <- sa[thisWav, ]
                            if(any(thisSa$Status == 'Start')) {
                                wavMap$startSample[w] <- 1
                            }
                        }
                    } else {
                        wavMap$startSample <- 1
                    }
                }
            }

            isMapped <- c(isMapped, folder[d])
        }
        dbMap[[d]] <- wavMap
    }

    # allFiles <- unique(sapply(dbMap, function(d) d$file))
    allFiles <- bind_rows(dbMap, .id = 'db')
    # combine with old then re-check for consecutive files within each db
    sameFile <- files(x)$recordings$file[files(x)$recordings$file %in% allFiles$file]
    for(f in seq_along(sameFile)) {
        whichNew <- allFiles$file == sameFile[f]
        newNA <- recNA(allFiles[whichNew, ])
        if(!newNA) {
            next
        }
        whichOrig <- files(x)$recordings$file == sameFile[f]
        origNA <- recNA(files(x)$recordings[whichOrig, ])
        if(!origNA) {
            allFiles[whichNew, ] <- files(x)$recordings[whichOrig, colnames(allFiles)]
        }
    }
    allFiles <- bind_rows(allFiles, files(x)$recordings)
    allFiles <- bind_rows(lapply(split(allFiles, allFiles$db), function(x) {
        checkConsecutive(x)
    }))
    allFiles <- bind_rows(lapply(split(allFiles, allFiles$fileGroup), function(x) {
        if(!is.na(x$startSample[1]) &&
           nrow(x) > 1 &&
           all(is.na(x$startSample[2:nrow(x)]))) {
            x$startSample[2:nrow(x)] <- cumsum(x$sampleLength[1:(nrow(x)-1)])
        }
        x
    }))
    if(any(is.na(allFiles$startSample))) {
        pamWarning('Unable to find "startSample" for some wav files, ',
                   'wav folder must contain all files processed if "Merge Contiguous Files" ',
                   'was selected in PAMGuard.')
    }
    files(x)$recordings <- allFiles
    # for(e in seq_along(events(x))) {
    #     thisFiles <- dbMap[[files(x[[e]])$db]]
    #     if(is.null(thisFiles) ||
    #        isFALSE(thisFiles)) next
    #     files(x[[e]])$recordings <- bind_rows(
    #         files(x[[e]])$recordings,
    #         thisFiles
    #     )
    #     files(x[[e]])$recordings <- checkConsecutive(files(x[[e]])$recordings)
    # }
    files(x)$recordings$numStart <- as.numeric(files(x)$recordings$start)
    files(x)$recordings$numEnd <- as.numeric(files(x)$recordings$end)
    x <- .addPamWarning(x)
    x
}
#' @export
#' @rdname addRecordings
#'
mapWavFolder <- function(folder=NULL, log=NULL, progress=TRUE) {
    if(is.null(folder)) {
        folder <- tk_choose.dir(caption = 'Select a folder containing your recording files.',
                                   default = getwd())
    }
    if(is.na(folder)) {
        return(NULL)
    }
    if(!dir.exists(folder)) {
        pamWarning('Provided folder ', folder, ' does not exist.')
        return(NULL)
    }
    folder <- normalizePath(folder, winslash = '/')
    wavs <- list.files(folder, full.names=TRUE, pattern = '\\.wav$', recursive=TRUE)
    if(length(wavs) == 0) {
        pamWarning('No wav files found in folder ', folder)
        return(NULL)
    }
    wavMap <- wavsToRanges(wavs, log, progress)
    if(is.null(wavMap) ||
       nrow(wavMap) == 0) {
        return(NULL)
    }
    wavMap <- checkConsecutive(wavMap)
    wavMap <- arrange(wavMap, .data$start)
    wavMap
}

wavsToRanges <- function(wav, log, progress=TRUE) {
    FOUNDFORMAT <- NULL
    if(progress) {
        cat('Mapping wav folder...\n')
        pb <- txtProgressBar(min=0, max=length(wav), style = 3)
        wix <- 1
    }
    noRead <- character(0)
    noReadMsg <- character(0)
    wavMap <- bind_rows(lapply(wav, function(x) {
        header <- try(readWave(x, header = TRUE), silent=TRUE)
        if(inherits(header, 'try-error')) {
            noRead <<- c(noRead, x)
            noReadMsg <<- c(noReadMsg, attr(header, 'condition')$message)
            len <- NA
            sampleLength <- NA
            sr <- NA
        } else {
            sr <- header$sample.rate
            len <- header$samples / sr
            sampleLength <- header$samples
        }
        format <- c(FOUNDFORMAT, c('pamguard', 'pampal', 'soundtrap', 'sm3', 'icListens1', 'icListens2'))
        for(f in format) {
            switch(
                f,
                'pamguard' = {
                    date <- gsub('.*([0-9]{8}_[0-9]{6}_[0-9]{3})\\.wav$', '\\1', x)
                    posix <- as.POSIXct(substr(date, 1, 15), tz = 'UTC', format = '%Y%m%d_%H%M%S')
                    if(is.na(posix)) next
                    millis <- as.numeric(substr(date, 17, 19)) / 1e3
                    if(!is.na(posix)) {
                        FOUNDFORMAT <<- f
                        break
                    }
                },
                'pampal' = {
                    date <- gsub('.*([0-9]{14}_[0-9]{3})\\.wav$', '\\1', x)
                    posix <- as.POSIXct(substr(date, 1, 14), tz = 'UTC', format = '%Y%m%d%H%M%S')
                    if(is.na(posix)) next
                    millis <- as.numeric(substr(date, 16, 18)) / 1e3
                    if(!is.na(posix)) {
                        FOUNDFORMAT <<- f
                        break
                    }
                },
                'soundtrap' = {
                    date <- gsub('.*\\.([0-9]{12})\\.wav$', '\\1', x)
                    posix <- as.POSIXct(date, format = '%y%m%d%H%M%S', tz='UTC')
                    millis <- 0
                    if(!is.na(posix)) {
                        FOUNDFORMAT <<- f
                        break
                    }
                },
                'sm3' = {
                    date <- gsub('.*\\_([0-9]{8}_[0-9]{6})Z?\\.wav$', '\\1', x)
                    posix <- as.POSIXct(date, format = '%Y%m%d_%H%M%S', tz='UTC')
                    millis <- 0
                    if(!is.na(posix)) {
                        FOUNDFORMAT <<- f
                        break
                    }
                },
                'icListens1' = {
                    date <- gsub('.*_([0-9]{8}-[0-9]{6})\\.wav$', '\\1', x)
                    posix <- as.POSIXct(date, format = '%Y%m%d-%H%M%S', tz='UTC')
                    millis <- 0
                    if(!is.na(posix)) {
                        FOUNDFORMAT <<- f
                        break
                    }
                },
                'icListens2' = {
                    date <- gsub('.*_([0-9]{6}-[0-9]{6})\\.wav$', '\\1', x)
                    posix <- as.POSIXct(date, format = '%y%m%d-%H%M%S', tz='UTC')
                    millis <- 0
                    if(!is.na(posix)) {
                        FOUNDFORMAT <<- f
                        break
                    }
                }
            ) #'AMAR668.9.20210823T231318Z.wav'
        }
        if(progress) {
            setTxtProgressBar(pb, value=wix)
            wix <<- wix + 1
        }
        # if(is.na(posix)) {
        #     pamWarning('Could not convert the name of the wav file ',
        #             x, ' to time properly.')
        #     return(NULL)
        # }
        if(!is.null(FOUNDFORMAT) &&
           FOUNDFORMAT == 'soundtrap') {
            hasLog <- which(gsub('\\.wav$', '', basename(x)) == log$file)
            if(length(hasLog) == 1) {
                len <- len + log$gap[hasLog]
                sampleLength <- sampleLength + log$gap[hasLog] * sr
            }
        }
        rng <- c(0, len) + posix + millis
        # if(any(is.na(rng))) {
        #     return(NULL)
        # }

        list(start=rng[1], end=rng[2], file=x, length=len, sampleLength = sampleLength, sr=sr)
    }))
    # if(length(badWav) > 0) {
    #     pamWarning('Unable to read wav files ', badWav, ' these are possibly corrupt.', n=6)
    # }
    noConvert <- is.na(wavMap$start)
    if(any(noConvert)) {
        pamWarning('Could not convert wav names to time properly for files ',
                   wavMap$file[noConvert], n=6)
    }
    if(length(noRead) > 0) {
        pamWarning('Could not read wav files ', noRead, n=6)
        pamWarning('Read failed with error messages ', noReadMsg, n=6)
    }
    wavMap
}
getSoundtrapLog <- function(log) {
    if(is.null(log) ||
       isFALSE(log) ||
       is.na(log) ||
       log == 'FALSE') {
        # log <- choose.dir(caption='Select a folder of Soundtrap log files (optional)')
        # cat('Select a folder of SoundTrap log files (optional)')
        # log <- tk_choose.dir(caption = 'Select a folder of SoundTrap log files (optional)',
        #                      default = getwd())
        return(data.frame(gap=0, sample=1, file='DNE'))
    }
    if(!dir.exists(log)) {
        pamWarning('Log folder ', log, ' does not exist.')
        return(data.frame(gap=0, sample=1, file='DNE'))
    }

    xFold <- list.files(log, full.names = TRUE, pattern = 'xml')
    if(length(xFold) == 0) {
        pamWarning('Log folder ', log, ' has no XML log files.')
        return(data.frame(gap=0, sample=1, file='DNE'))
    }
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
        list(micros=0, sample=1, file = gsub('\\.log\\.xml$', '', basename(x)))
    })
    missing <- bind_rows(missing) %>%
        group_by(.data$file) %>%
        summarise(gap = sum(.data$micros)/1e6) %>%
        ungroup()
    missing
}

# check for consecutive recordings and mark with single number
checkConsecutive <- function(wavMap) {
    wavMap <- distinct(wavMap, .data$file, .keep_all=TRUE)
    wavMap <- arrange(wavMap, .data$start)
    wavMap$fileGroup <- 1
    wavMap$timeDiff <- 0
    if(nrow(wavMap) > 1) {
        fg <- 1
        for(i in 2:nrow(wavMap)) {
            wavMap$timeDiff[i] <- as.numeric(difftime(wavMap$start[i], wavMap$end[i-1], units='secs'))
            # if theres a tiny negative time difference because ST files are weird
            # then just make them equal because they should be
            # if(wavMap$timeDiff[i] < 0 &&
            #    wavMap$timeDiff[i] > -.05) {
            #     # wavMap$end[i-1] <- wavMap$start[i]
            #     wavMap$timeDiff[i] <- 0
            # }
            # if(wavMap$timeDiff[i] != 0) {
            # # if(wavMap$end[i-1] != wavMap$start[i]) {
            #     fg <- fg + 1
            # }
            if(is.na(wavMap$timeDiff[i]) ||
               wavMap$timeDiff[i] > .02 ||
               wavMap$timeDiff[i] < -.05) {
                fg <- fg + 1
            }
            wavMap$fileGroup[i] <- fg
        }
    }
    wavMap
}
