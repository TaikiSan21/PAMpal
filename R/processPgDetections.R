#' @title Load and Process Detections from Pamguard
#'
#' @description Loads and processes acoustic detection data that has been
#'   run through Pamguard. Uses the binary files and database(s) contained
#'   in \code{pps}, and will group your data into events by the
#'   grouping present in the 'OfflineEvents' and 'Detection Group Localiser'
#'   tables (\code{mode = 'db'}) or by the grouping specified by start and end
#'   times in the supplied \code{grouping} (\code{mode = 'time'}), or by start and
#'   end of recording files (\code{mode = 'recording'}). Will apply
#'   all processing functions in \code{pps} to the appropriate modules
#'
#' @param pps a \linkS4class{PAMpalSettings} object containing the databases,
#'   binaries, and functions to use for processing data. See
#'   \code{\link[PAMpal]{PAMpalSettings}}. Can also be an \linkS4class{AcousticStudy}
#'   object, in which case the \code{pps} slot will be used.
#' @param mode selector for how to organize your data in to events. \code{db}
#'   will organize by events based on tables in the databases. \code{time}
#'   will organize into events based on timestamps provided in \code{grouping}.
#'   \code{recording} will organize events by the start and end times of recording
#'   files found in the database. For \code{time} and \code{recording}, ALL detections
#'   between the start and end times are included, for \code{db} only selected
#'   detections are included.
#' @param id an event name or id for this study, will default to today's date if
#'   not supplied (recommended to supply your own informative id)
#' @param grouping For \code{mode = 'db'}, the table to group events by.
#'   Either \code{event} to use the OfflineEvents table, or \code{detGroup} to
#'   use the detection group localiser module groups.
#'
#'   For \code{mode = 'time'},
#'   this should be a data frame with three mandatory columns and 1 row
#'   for each separate event. The mandatory columns are \code{start}, \code{end},
#'   and \code{id}. \code{start} and \code{end} should specify the
#'   start and end time of the event and must be in UTC. \code{id} should
#'   specify a unique id for each event. There are also optional columns
#'   \code{species}, \code{db}, and \code{sr}. \code{species} should provide a
#'   species ID if it is available. \code{db} and \code{sr} are the corresponding
#'   database and sample rate to associate with a particular event, these typically
#'   do not need to be specified as the function will attempt to automatically match
#'   them based on the times of the events and the databases. Note that \code{db}
#'   must be the full filepath to the database. If a clear match is not found then
#'   the user will be prompted to either select from a list or input the proper
#'   sample rate.
#'
#'   \code{grouping} can be supplied either as a data frame or as
#'   a filepath to a csv file.
#' @param format the date format for the \code{start} and \code{end} columns
#'   in \code{grouping} if it is a csv. Times are assumed to be UTC. See
#'   details section of \link{strptime} for more information on how to properly
#'   format this
#' @param progress logical flog to show progress bars
#' @param verbose logical flag to show messages
#' @param \dots additional arguments to pass onto to different methods
#'
#' @details If \code{mode} is not specified, it will try to be automatically determined
#'   in the following order. 1) if a \code{grouping} data.frame or CSV is provided, then
#'   \code{mode='time'} will be used. 2) If there are labelled events present in the
#'   database, \code{mode='db'} will be used. 3) \code{mode='recording'} will be used,
#'   which should be equivalent to loading all possible data.
#'
#' @return an \linkS4class{AcousticStudy} object with one \linkS4class{AcousticEvent}
#'   for each event in the \code{events} slot, and the \linkS4class{PAMpalSettings} object
#'   used stored in the \code{pps} slot.
#'
#' @author Taiki Sakai \email{taiki.sakai@@noaa.gov}
#'
#' @examples
#'
#' exPps <- new('PAMpalSettings')
#' exPps <- addDatabase(exPps, system.file('extdata', 'Example.sqlite3', package='PAMpal'))
#' exPps <- addBinaries(exPps, system.file('extdata', 'Binaries', package='PAMpal'))
#' exClick <- function(data) {
#'     standardClickCalcs(data, calibration=NULL, filterfrom_khz = 0)
#' }
#' exPps <- addFunction(exPps, exClick, module = 'ClickDetector')
#' exPps <- addFunction(exPps, roccaWhistleCalcs, module='WhistlesMoans')
#' exPps <- addFunction(exPps, standardCepstrumCalcs, module = 'Cepstrum')
#' # process events labelled within the Pamguard database
#' exStudyDb <- processPgDetections(exPps, mode='db', id='Example')
#' # can also give an AcousticStudy as input and it will use same functions and data
#' reprocess <- processPgDetections(exStudyDb, mode='db', id='Reprocess')
#' # process events with manually set start/end times
#' grp <- data.frame(start = as.POSIXct('2018-03-20 15:25:10', tz='UTC'),
#'                   end = as.POSIXct('2018-03-20 15:25:11', tz='UTC'),
#'                   id = 'GroupExample')
#' exStudyTime <- processPgDetections(exPps, mode='time', grouping=grp, id='Time')
#' # process events by recording event
#' exStudyRecording <- processPgDetections(exPps, mode='recording', id='Recording')
#'
#' @importFrom PamBinaries loadPamguardBinaryFile
#' @importFrom PAMmisc squishList
#' @importFrom RSQLite dbConnect dbListTables dbReadTable dbDisconnect SQLite dbListFields
#' @importFrom stringr str_trim
#' @importFrom tcltk tk_choose.files
#' @importFrom purrr transpose
#' @import dplyr
#' @export
#'
processPgDetections <- function(pps, mode = c('db', 'time', 'recording'), id=NULL, grouping=NULL,
                                format='%Y-%m-%d %H:%M:%OS', progress=TRUE, verbose=TRUE, ...) {
    # auto check for mode
    if(missing(mode)) {
        mode <- autoMode(pps, grouping)
        if(verbose) {
            cat('Processing mode not provided, running with mode="',
                mode, '"\n', sep='')
        }
    }
    mode <- match.arg(mode)
    if(is.AcousticStudy(pps)) {
        if(mode == 'time' &&
           is.null(grouping) &&
           !is.null(ancillary(pps)$grouping)) {
            if(verbose) {
                cat('Found a grouping file in the provided AcousticStudy object,',
                    'to use a different grouping file specify with the grouping argument.')
            }
            grouping <- ancillary(pps)$grouping
        }
        pps <- pps(pps)
    }
    if(!is.PAMpalSettings(pps)) {
        stop(paste0(pps, ' is not a PAMpalSettings object. Please create one with',
                    ' function "PAMpalSettings()"'), call.=FALSE)
    }
    result <- switch(mode,
                     'db' = processPgDb(pps=pps, grouping=grouping, id=id,
                                                  progress=progress, ...),
                     'time' = processPgTime(pps=pps, grouping=grouping, format=format, id=id,
                                                      progress=progress),
                     'recording' = {
                         grouping <- bind_rows(lapply(pps@db, function(x) {
                             wavToGroup(x)
                         }))
                         if(is.null(grouping) ||
                            nrow(grouping) == 0) {
                             warning('Unable to load Sound_Acquisition data properly, cannot run mode="recording"', call.=FALSE)
                             return(NULL)
                         }
                         processPgTime(pps=pps, grouping=grouping, format=format, id=id,
                                                 progress=progress)
                     }
    )
    checkStudy(result)
    result <- .addPamWarning(result)
    result
}

# ---- separate methods ----

#' @importFrom utils setTxtProgressBar txtProgressBar
#'
processPgTime <- function(pps, grouping=NULL, format='%Y-%m-%d %H:%M:%OS', id=NULL,
                                    progress=progress) {
    # start with checking grouping - parse csv if missing or provided as character and fmt times
    grouping <- checkGrouping(grouping, format)
    # this is a flag to see if any manual entries happened to grouping
    editGroup <- FALSE
    binList <- pps@binaries$list
    binFuns <- pps@functions
    allDbs <- pps@db

    # study <- AcousticStudy(id=id, pps = pps)
    # Check for what DB shit should be associated with, get full list of SA data
    # first, gonna match event times to that since its roughly the times assoicated
    # with a database
    saList <- lapply(allDbs, readSa)
    names(saList) <- allDbs

    if(!('db' %in% colnames(grouping))) {
        grouping$db <- NA_character_
    }
    # if they are there and are valid, assume they assigned
    dbToAssign <- which(sapply(basename(grouping$db), function(d) {
        if(is.na(d)) {
            return(TRUE)
        }
        !any(grepl(d, allDbs))
    }))
    # match db to events
    if(length(allDbs) == 1) {
        grouping$db[dbToAssign] <- allDbs
    } else {
        for(i in dbToAssign) {
            if(is.na(grouping$db[i]) ||
               !any(grepl(grouping$db[i], allDbs, fixed=TRUE))) {
                dbPossible <- allDbs[sapply(saList, function(x) {
                    inInterval(c(grouping$start[i], grouping$end[i]), x)
                })]
            } else { # case if you just specified basename of the database it will find it
                dbPossible <- grep(grouping$db[i], allDbs, value=TRUE, fixed=TRUE)
            }

            if(length(dbPossible) == 0 ||
               is.na(dbPossible)) {
                editGroup <- TRUE
                myTitle <- paste0('No matching database found for event ', grouping$id[i],
                                  ' based on times, please choose one or select "0" to',
                                  ' leave as NA.')
                myChoice <- menu(title = myTitle, choices = c(allDbs, 'Exit function call (no processing will occur)'))
                if(myChoice == length(allDbs) + 1) {
                    stop('Exiting function call', call.=FALSE)
                }
                if(myChoice == 0) {
                    dbPossible <- NA_character_
                } else {
                    dbPossible <- allDbs[myChoice]
                }
            } else if(length(dbPossible) > 1) {
                editGroup <- TRUE
                myTitle <- paste0('Multiple candidate datbases found for event "', grouping$id[i],
                                  '" based on times, select one to associate with this event.')
                myChoice <- menu(title = myTitle, choices = dbPossible)
                if(myChoice == 0) {
                    dbPossible <- NA_character_
                } else {
                    dbPossible <- dbPossible[myChoice]
                }
            }
            grouping$db[i] <- dbPossible
        }
    }
    failBin <- 'Havent started'
    on.exit({
        # only do the saving if anything had to be done by the user
        if(editGroup) {
            time <- gsub(' ', '_', as.character(Sys.time()))
            time <-gsub(':', '-', time)
            fileName <- paste0(time, '_GroupingData.Rdata')
            message('\nOops! It looks like something went wrong and the function ',
                'stopped before finishing. Your "grouping"',
                ' data has been saved in the current working directory as:\n',
                '   ', fileName, '\nYou can supply this to "grouping" next time you ',
                'run getPgDetections to avoid re-selecting options with:\n',
                '   newGrouping <- readRDS("', fileName, '")')
            saveRDS(grouping, file = fileName)
        }
        message('\nLast file I tried to read: ', failBin)
    })

    if(!('sr' %in% colnames(grouping))) {
        grouping$sr <- NA_integer_
    }

    # assign each db in grouping to its unique SRs so we dont have to search again later
    # saByDb <- lapply(saList, function(x) unique(x$sampleRate))
    # for(d in 1:nrow(grouping)) {
    #     if(is.na(grouping$db[d])) next
    #     grouping$sr[d] <- saByDb[grouping$db[d]]
    # }
    grouping <- bind_rows(lapply(split(grouping, grouping$db), function(d) {
        thisSa <- saList[[unique(d$db)]]
        if(is.null(thisSa) ||
           nrow(thisSa) == 0) {
            return(d)
        }
        d$UTC <- d$end
        d <- matchSR(d, thisSa)
        d$UTC <- NULL
        d$sr <- d$sampleRate
        d$sampleRate <- NULL
        d
    }))
    # were gonna match SR by database, only need manual input if we have any
    # missing DBs
    if(any(is.na(grouping$sr))) {
        editGroup <- TRUE
        cat('\nCould not read sample rate from database,',
            'what sample rate should be used for these events?\n')
        sr <- readline()
        grouping$sr[is.na(grouping$sr)] <- as.numeric(sr)
    }

    # from here can check "simple SR" mode - all SR in DBs and
    # the one we selected are the same, avoid doing shit later

    calibrationUsed <- names(pps@calibration[[1]])
    if(length(calibrationUsed)==0) calibrationUsed <- 'None'

    binExists <- file.exists(binList)
    if(sum(binExists) == 0) {
        stop('No valid binary files found. Either none have been added, or the ',
             'path has changed or is incorrect. Please add again with function ',
             '"addBinaries".', call.=FALSE)
    }
    if(any(!binExists)) {
        contChoice <- menu(title=paste0(sum(!binExists), ' out of ', length(binExists),
                                        ' binary files could not be found, would you',
                                        ' like to continue processing or stop to investigate?'),
                           choices=c('Continue', 'Stop'))
        if(contChoice == 2) {
            stop('Stopping, no processing has been done', call.=FALSE)
        }
    }
    binList <- binList[binExists]
    if(progress) {
        cat('Processing binary files... \n')
        pb <- txtProgressBar(min=0, max=length(binList), style=3)
    }
    modList <- c('ClickDetector', 'WhistlesMoans', 'Cepstrum', 'GPLDetector')
    modWarn <- c(FALSE, FALSE, FALSE, FALSE)
    names(modWarn) <- modList
    binData <- lapply(binList, function(bin) {
        # should i do here - read in head/foot only, then check those
        # times against grouplist, if none can skip, if one we know
        # what db to match sr with. if more than one... hope they have the
        # same SR? or go fys?

        # debugger
        failBin <<- bin
        # flag if weve loaded data, need because incomplete binaries dont have footer for check
        loaded <- FALSE
        thisHFOnly <- loadPamguardBinaryFile(bin, skipData=TRUE)$fileInfo
        # if either of these isnt present we need to load binary file completely so check first
        dateBounds <- c(thisHFOnly$fileHeader$dataDate, thisHFOnly$fileFooter$dataDate)
        if(length(dateBounds) == 2) {
            binBounds <- convertPgDate(dateBounds)
        } else {
            thisBin <- loadPamguardBinaryFile(bin)
            loaded <- TRUE
            dataLen <- length(thisBin$data)
            if(dataLen == 0) {
                if(progress) {
                    setTxtProgressBar(pb, value=which(binList==bin))
                }
                return(NULL)
            }
            binBounds <- convertPgDate(c(thisBin$data[[1]]$date, thisBin$data[[dataLen]]$date))
        }

        evPossible <- (binBounds[1] >= grouping$start & binBounds[1] <= grouping$end) |
            (binBounds[2] >= grouping$start & binBounds[2] <= grouping$end) |
            (binBounds[1] <= grouping$start & binBounds[2] >= grouping$end)

        # if not overlapping any events, skip doing data part mobetta
        if(!any(evPossible)) {
            if(progress) {
                setTxtProgressBar(pb, value=which(binList==bin))
            }
            return(NULL)
        }
        if(!loaded) {
            thisBin <- loadPamguardBinaryFile(bin)
        }
        if(length(thisBin$data) == 0) {
            if(progress) {
                setTxtProgressBar(pb, value=which(binList==bin))
            }
            return(NULL)
        }
        modType <- getModuleType(thisBin)
        if(length(binFuns[[modType]]) == 0) {
            if(isFALSE(modWarn[[modType]])) {
                modWarn[[modType]] <<- TRUE
                pamWarning('No functions for processing Module Type: ', modType)
            }
        }
        srPossible <- unique(unlist(grouping$sr[evPossible]))
        if(length(srPossible) == 1) {
            for(i in seq_along(thisBin$data)) {
                thisBin$data[[i]]$sr <- srPossible
            }
        } else if(length(srPossible) > 1) {
            evDbs <- unique(grouping$db[evPossible])
            thisSa <- do.call(rbind, saList[evDbs])
            binTimes <- bind_rows(lapply(thisBin$data, function(x) {
                list(UID = x$UID, UTC = x$date)
            }))
            binTimes$UTC <- convertPgDate(binTimes$UTC)
            binTimes <- matchSR(binTimes, thisSa)
            for(i in seq_along(thisBin$data)) {
                thisBin$data[[i]]$sr <- binTimes$sampleRate[i]
            }
        }
        thisBinData <- calculateModuleData(thisBin, binFuns, pps@settings)
        if(progress) {
            setTxtProgressBar(pb, value=which(binList==bin))
        }
        thisBinData
    })
    if(progress) {
        cat('\n') # space after progress bar finished
    }
    binData <- binData[sapply(binData, function(x) !is.null(x))]
    if(length(binData) == 0) {
        stop(paste0('None of the binary files contained data for any of the events.',
                    ' Please check that times are in UTC and the correct binary folder was supplied.'), call.=FALSE)
    }
    # for clicks we have split the broad detector into separate ones by classification
    binData <- lapply(binData, function(x) split(x, x$detectorName))
    binData <- unlist(binData, recursive = FALSE)
    binData <- squishList(binData)

    acousticEvents <- vector('list', length = nrow(grouping))
    evName <- as.character(grouping$id)

    colsToDrop <- c('Id', 'comment', 'sampleRate', 'detectorName', 'parentID',
                    'sr', 'callType', 'newUID')
    names(acousticEvents) <- evName
    noDetEvent <- character(0)
    for(i in seq_along(acousticEvents)) {
        thisData <- lapply(binData, function(x) {
            data <- filter(x, x$UTC >= grouping$start[i], x$UTC <= grouping$end[i])
            if(nrow(data) == 0) return(NULL)
            data
        })
        # Check possible DBs by start/end time of events in sa list earlier
        thisData <- thisData[sapply(thisData, function(x) !is.null(x))]
        binariesUsed <- sapply(thisData, function(x) unique(x$BinaryFile)) %>%
            unlist(recursive = FALSE) %>% unique()
        # binariesUsed <- unlist(sapply(binariesUsed, function(x) grep(x, binList, value=TRUE), USE.NAMES = FALSE))
        # Check and warning here for empty event
        if(length(thisData) == 0) {
            noDetEvent <- c(noDetEvent, names(acousticEvents)[i])
        }
        thisData <- lapply(thisData, function(x) {
            x$BinaryFile <- basename(x$BinaryFile)
            thisType <- unique(x$callType)
            x <- dropCols(x, colsToDrop)
            attr(x, 'calltype') <- thisType
            x
        })
        thisSr <- grouping$sr[[i]]
        if(is.na(grouping$db[i])) {
            thisSource <- 'Not Found'
        } else {
            filtSa <- saList[[grouping$db[i]]]
            if(is.null(filtSa)) {
                thisSource <- 'Not Found'
            } else {
                filtSa <- filter(filtSa, filtSa$UTC <= grouping$end[i], filtSa$UTC >= grouping$start[i])
                thisSource <- unique(filtSa$SystemType)
            }
        }
        thisDb <- allDbs[grepl(basename(grouping$db[i]), allDbs)]
        acousticEvents[[i]] <-
            AcousticEvent(id=evName[i], detectors = thisData, settings = list(sr = thisSr, source = thisSource),
                          files = list(binaries=binariesUsed, db=thisDb, calibration=calibrationUsed))
    }
    if(length(noDetEvent) > 0) {
        pamWarning('No detections in Event(s) ', noDetEvent, n=6)
    }
    if('species' %in% colnames(grouping)) {
        grouping$species <- as.character(grouping$species)
        acousticEvents <- setSpecies(acousticEvents, method = 'manual', value = grouping$species)
    }
    allDbs <- unique(unlist(lapply(acousticEvents, function(x) {
        files(x)$db
    })))
    allBins <- unique(unlist(lapply(acousticEvents, function(x) {
        files(x)$binaries
    })))
    study <- AcousticStudy(id=id, events = acousticEvents, pps = pps,
                           files = list(db=allDbs, binaries=allBins))
    settings(study)$detectors <- pps@settings$detectors
    # events(study) <- acousticEvents
    # files(study) <- list(db=allDbs, binaries=allBins)
    ancillary(study)$grouping <- grouping
    # study <- .addPamWarning(study)
    on.exit() # this cancels the on.exit 'save my grouping' call that is there if you crash
    study
}

#'
processPgDb <- function(pps, grouping=c('event', 'detGroup'), id=NULL,
                                  progress=TRUE, ...) {
    allDb <- pps@db

    # awk diff init values between modes have to reset this here
    if(is.null(grouping)) {
        grouping <- c('event', 'detGroup')
    }

    nBin <- sum(sapply(allDb, nBins))
    if(nBin == 0) {
        warning('No detections found within database, are you sure you want',
                'to run with mode="db" ?')
        return(NULL)
    }
    if(progress) {
        cat('Processing databases... \n')
        pb <- txtProgressBar(min=0, max=nBin, style=3)
    }
    binNo <- 1
    modList <- c('ClickDetector', 'WhistlesMoans', 'Cepstrum', 'GPLDetector')
    modWarn <- c(FALSE, FALSE, FALSE, FALSE)
    tarMoCols <- ppVars()$tarMoCols
    names(modWarn) <- modList
    allAcEv <- lapply(allDb, function(db) {
        tryCatch({
            missBins <- character(0)
            failBin <- 'No file processed'
            binList <- pps@binaries$list
            binFuns <- pps@functions
            dbData <- getDbData(db=db, grouping=grouping, extraCols=tarMoCols, ...)
            if(is.null(dbData) ||
               nrow(dbData) == 0) {
                pamWarning('No detections found in database ',
                        basename(db), '.')
                # setTxtProgressBar(pb, value = evNo)
                # evNo <- evNo + 1
                return(NULL)
            }
            thisSr <- unique(dbData$sampleRate)
            if(length(thisSr) > 1) {
                pamWarning('More than 1 sample rate found in database ',
                        basename(db),'.')
            }
            thisSource <- unique(dbData$SystemType)
            dbData <- select(dbData, -.data$SystemType)
            calibrationUsed <- names(pps@calibration[[1]])
            if(length(calibrationUsed)==0) calibrationUsed <- 'None'

            dbData <- lapply(
                split(dbData, dbData$BinaryFile), function(x) {
                    if(progress) {
                        setTxtProgressBar(pb, value = binNo)
                    }
                    binNo <<- binNo + 1
                    failBin <<- x$BinaryFile[1]
                    thisBin <- getMatchingBinaryData(x, binList, basename(db))
                    if(length(thisBin)==0) {
                        missBins <<- c(missBins, x$BinaryFile[1])
                        # warning('Could not find the matching binary file for ', x$BinaryFile[1],
                        #         ' in database ', basename(db), call.=FALSE)
                        return(NULL)
                    }
                    modType <- getModuleType(thisBin)
                    if(length(binFuns[[modType]]) == 0) {
                        if(isFALSE(modWarn[[modType]])) {
                            modWarn[[modType]] <<- TRUE
                            pamWarning(paste0('No functions for processing Module Type: ', modType))
                        }
                    }
                    binData <- calculateModuleData(thisBin, binFuns, pps@settings)
                    if(!is.null(binData)) {
                        noMatch <- which(!(x$UID %in% binData$UID))
                        x$UID[noMatch] <- x$newUID[noMatch]
                        x$BinaryFile <- NULL
                        binData %>%
                            # select(-.data$BinaryFile) %>%
                            inner_join(x, by='UID') %>%
                            distinct()
                    }
                }
            )
            if(length(missBins) > 0) {
                pamWarning('Could not find the matching binary files for binaries ',
                                  missBins,
                                  ' in database ', basename(db), n=3)
            }
            # This is a list for each binary, we want for each detector
            dbData <- dbData[sapply(dbData, function(x) !is.null(x))]

            dbData <- lapply(dbData, function(x) split(x, x$detectorName))
            names(dbData) <- NULL
            dbData <- unlist(dbData, recursive = FALSE)
            dbData <- squishList(dbData)

            # Split into events, then swap from Detector(Events) to Event(Detectors)
            # .names necessary to make sure we have all event numbers
            dbData <- transpose(
                lapply(dbData, function(x) split(x, x$parentID)),
                .names = unique(unlist(sapply(dbData, function(x) x$parentID)))
            )

            # Should this function store the event ID? Right now its just the name
            # in the list, but is this reliable? Probably not

            colsToDrop <- c('Id', 'comment', 'sampleRate', 'detectorName', 'parentID',
                            'sr', 'callType', 'newUID', tarMoCols)

            acousticEvents <- lapply(dbData, function(ev) {
                ev <- ev[sapply(ev, function(x) !is.null(x))]
                binariesUsed <- sapply(ev, function(x) unique(x$BinaryFile)) %>%
                    unlist(recursive = FALSE) %>% unique()
                binariesUsed <- unlist(sapply(binariesUsed, function(x) grep(x, binList, value=TRUE, fixed=TRUE), USE.NAMES = FALSE))

                evId <- paste0(gsub('\\.sqlite3', '', basename(db)), '.', unique(ev[[1]]$parentID))
                if(all(tarMoCols %in% colnames(ev[[1]]))) {
                    evTarMo <- ev[[1]][1, tarMoCols]
                } else {
                    evTarMo <- data.frame(matrix(NA, nrow=1, ncol=length(tarMoCols)))
                    colnames(evTarMo) <- tarMoCols
                }

                evComment <- unique(ev[[1]]$comment)
                if(is.null(evComment)) {
                    evComment <- NA
                }
                ev <- lapply(ev, function(x) {
                    x$BinaryFile <- basename(x$BinaryFile)
                    thisType <- unique(x$callType)
                    x <- dropCols(x, colsToDrop)
                    attr(x, 'calltype') <- thisType
                    x
                })
                acEv <- AcousticEvent(id = evId, detectors = ev, settings = list(sr = thisSr, source=thisSource),
                              files = list(binaries=binariesUsed, db=db, calibration=calibrationUsed),
                              localizations = list(PGTargetMotion = evTarMo))
                ancillary(acEv)$eventComment <- evComment
                acEv
            })
            # setTxtProgressBar(pb, value = evNo)
            # evNo <- evNo + 1
            acousticEvents
        },
        error = function(e) {
            message('\nError in processing db ', basename(db), ' during binary file ', failBin)
            message('\nError message:\n')
            print(e)
            # setTxtProgressBar(pb, value = evNo)
            # evNo <- evNo + 1
            return(NULL)
        })
    })

    if(progress) {
        cat('\n')
    }
    names(allAcEv) <- gsub('\\.sqlite3', '', basename(allDb))
    allAcEv <- unlist(allAcEv, recursive = FALSE)
    allDbs <- unique(unlist(lapply(allAcEv, function(x) {
        files(x)$db
    })))
    allBins <- unique(unlist(lapply(allAcEv, function(x) {
        files(x)$binaries
    })))
    on.exit()
    study <- AcousticStudy(id=id, events = allAcEv, pps = pps,
                  files = list(db=allDbs, binaries=allBins))
    # study <- .addPamWarning(study)
    settings(study)$detectors <- pps@settings$detectors
    study
}

# ---- not exported helpers ----
getDbData <- function(db, grouping=c('event', 'detGroup'), label=NULL, extraCols = NULL, doSR=TRUE) {
    # Combine all click/event tables, even by diff detector. Binary will have det name
    con <- dbConnect(SQLite(), db)
    on.exit(dbDisconnect(con))
    tables <- dbListTables(con)
    # Read in event data from either offlineclicks/events or detection
    # group localiser. Click version has common naming convention,
    # det group does not so we have to go look it up. If we are just
    # reading in all the data we only care about SA data
    if(is.null(grouping)) {
        grouping <- c('event', 'detGroup')
    }
    if(length(grouping) > 1) {
        return(
            suppressPamWarnings(
                bind_rows(
                    lapply(grouping, function(x) {
                        getDbData(db, x, label, extraCols, doSR)
                    }))
            )
        )
    }
    switch(match.arg(grouping),
           'event' = {
               detTables <- grep('OfflineClicks', tables, value=TRUE)
               eventTables <- grep('OfflineEvents', tables, value=TRUE)
               # eventColumns <- c('UID', 'eventType', 'comment')
               if(is.null(label)) {
                   label <- 'eventType'
               }
               eventColumns <- c('Id', label, 'comment')
               evName <- 'OE'
           },
           'detGroup' = {
               modules <- dbReadTable(con, 'PamguardModules')
               dgTables <- modules %>%
                   mutate(Module_Name=str_trim(.data$Module_Name),
                          Module_Type=str_trim(.data$Module_Type)) %>%
                   filter(.data$Module_Name == 'Detection Group Localiser') %>%
                   distinct(.data$Module_Type, .data$Module_Name)
               dgNames <- gsub(' ',  '_', dgTables$Module_Type)
               detTables <- sapply(dgNames, function(x) grep(x, tables, value=TRUE))
               eventTables <- detTables[!grepl('Children', detTables)]
               detTables <- detTables[grepl('Children', detTables)]
               # eventColumns <- c('UID', 'Text_Annotation')
               if(is.character(eventTables)) {
                   dglCols <- dbListFields(con, eventTables[1])
                   label <- parseDglLabel(label, dglCols)
               } else {
                   label <- NA
               }

               eventColumns <- c('Id', label)
               evName <- 'DGL'
           },
           {
               stop("I don't know how to group by ", grouping, '.\n', call.=FALSE)
           }
    )

    if(length(detTables)==0 ||
       length(eventTables)==0) {
        pamWarning('Could not find event tables for grouping method "', grouping,
                '" in database ', basename(db), which = -1)
        return(NULL)
    }
    allDetections <- bind_rows(
        lapply(detTables, function(table) {
            dt <- dbReadTable(con, table)
            if(is.null(dt) ||
               nrow(dt) == 0) {
                return(NULL)
            }
            dt
        })
    )
    if(nrow(allDetections)==0) {
        pamWarning('No detections found for grouping method "', grouping,
                '" in database ', basename(db), which = -1)
        return(NULL)
    }

    allEvents <- bind_rows(
        lapply(eventTables, function(table) {
            et <- dbReadTable(con, table)
            if(is.null(et) ||
               nrow(et) == 0) {
                return(NULL)
            }
            et
        })
    )
    if(nrow(allEvents)==0) {
        pamWarning('No events found for grouping method "', grouping,
                '" in database ', basename(db), which = -1)
        return(NULL)
    }

    eventColumns <- eventColumns[eventColumns %in% colnames(allEvents)]
    allEvents <- select(allEvents, any_of(c(eventColumns, extraCols)))

    # Do i want all detections in clicks, or only all in events?
    # left_join all det, inner_join ev only
    if(!('Id' %in% names(allEvents)) ||
       !('parentID' %in% names(allDetections))) {
        message('Id and parentID columns not found in database ', basename(db),
                ', these are required to process data.')
        return(NULL)
    }

    allDetections <- inner_join(
        allDetections, allEvents, by=c('parentID'='Id')
    )
    if(!('newUID' %in% colnames(allDetections))) {
        allDetections$newUID <- -1
    }
    allDetections <- allDetections %>%
        mutate(BinaryFile = str_trim(.data$BinaryFile),
               # UTC = as.POSIXct(as.character(UTC), format='%Y-%m-%d %H:%M:%OS', tz='UTC')) %>%
               UTC = pgDateToPosix(.data$UTC)) %>%

        select(any_of(unique(c(eventColumns, 'UTC', 'UID', 'parentID', 'BinaryFile', 'newUID', extraCols))))

    # rename column to use as label - standardize across event group types
    colnames(allDetections)[which(colnames(allDetections)==label)] <- 'eventLabel'
    if(!('eventLabel' %in% colnames(allDetections))) {
        allDetections$eventLabel <- 'NOLABELFOUND'
    }
    if(doSR) {
        allDetections <- matchSR(allDetections, db, extraCols=c('SystemType'))
    }

    # apply str_trim to all character columns
    whichChar <- which(sapply(allDetections, function(x) 'character' %in% class(x)))
    for(i in whichChar) {
        allDetections[, i] <- str_trim(allDetections[, i])
    }
    allDetections <- select(allDetections, -.data$UTC)
    allDetections$UID <- as.character(allDetections$UID)
    allDetections$newUID <- as.character(allDetections$newUID)
    allDetections$parentID <- paste0(evName, allDetections$parentID)
    allDetections
}

getMatchingBinaryData <- function(dbData, binList, dbName, idCol = 'UID') {
    # dbData here has a single BinaryFile in it, we've split by that before here
    dbData <- arrange(dbData, .data[[idCol]])
    # This breaks if 'dbData' doesnt have binaryfile...
    # Borked if UID mismatch between dems
    binFile <- dbData$BinaryFile[1]
    allBinFiles <- grep(binFile, binList, value=TRUE, fixed=TRUE)
    if(length(allBinFiles)==0) {
        return(NULL)
    }
    dbData$matched <- FALSE
    for(bin in seq_along(allBinFiles)) {
        thisBin <- loadPamguardBinaryFile(allBinFiles[bin], keepUIDs = unique(c(dbData[['UID']], dbData[['newUID']])))

        # We've found the right file if theres any data
        if(length(thisBin$data) > 0) {
            # thisBin$data <- thisBin$data[names(thisBin$data) %in% dbData[[idCol]]]
            dbData$matched[dbData[['UID']] %in% names(thisBin$data)] <- TRUE
            dbData$matched[dbData[['newUID']] %in% names(thisBin$data)] <- TRUE
        }
        if(bin == 1) {
            result <- thisBin
        } else {
            result$data <- c(result$data, thisBin$data)
        }
        if(all(dbData$matched)) {
            break
        }
    }
    if(length(result$data) == 0) {
        return(NULL)
    }
    if(!all(dbData$matched)) {
        pamWarning('UID(s) ', dbData$UID[!dbData$matched],
                       ' are in databases ', dbName, ' but not in binary file ', binFile)
    }
    # do SR match after, if onyl one we can just simple match
    if(length(unique(dbData$sampleRate)) == 1) {
        for(i in seq_along(result$data)) {
            result$data[[i]]$sr <- dbData$sampleRate[1]
        }
    } else {
        for(i in names(result$data)) {
            ix <- min(which(dbData[['UID']] == i | dbData[['newUID']] == i))
            result$data[[i]]$sr <- dbData$sampleRate[ix]
        }
    }
    result
}

checkGrouping <- function(grouping, format) {
    if(is.null(grouping)) {
        # cat('Please provide a csv file with columns "start", "end", "id", and',
        #     'optionally "species" to group detections into events.')
        # grouping <- tk_choose.files(caption = 'Select event time csv file:', multi = FALSE)
        stop('Grouping file with columns "start", "end", and "id" must be provided to',
             ' group detections into events.')
    }
    if(inherits(grouping, 'character')) {
        if(!file.exists(grouping)) {
            stop('Provided grouping file does not exist, please provide a csv file with',
                ' columns "start", "end", and "id" to group detections into events.')
            # grouping <- tk_choose.files(caption = 'Select event time csv file:', multi = FALSE)
        }
        # grouping <- read_csv(grouping, col_types = cols(.default=col_character()))
        grouping <- read.csv(grouping, stringsAsFactors = FALSE, colClasses = 'character')
        if('sr' %in% colnames(grouping)) {
            grouping$sr <- as.numeric(grouping$sr)
        }
    }
    colsNeeded <- c('start', 'end', 'id')
    if(inherits(grouping, 'data.frame')) {
        colnames(grouping) <- tolower(colnames(grouping))
        if(!all(colsNeeded %in% colnames(grouping))) {
            stop('"grouping" must have columns "start", "end" and "id".')
        }
        # if times arent posix, convert and check that it worked
        if(!inherits(grouping$start, 'POSIXct') ||
           !inherits(grouping$end, 'POSIXct')) {

            grouping$start <- parseUTC(grouping$start, format)
            grouping$end <- parseUTC(grouping$end, format)

            if(any(is.na(grouping$start)) ||
               any(is.na(grouping$end))) {
                stop('Some event start/end times were not able to be converted, please check to',
                        ' ensure that format argument matches date format in provided file.')
            }
            # checkDate <- menu(title = paste0('\nThe first event start time is ', grouping$start[1],
            #                                  ', does this look okay?'),
            #                   choices = c('Yes, continue processing.',
            #                               "No. I'll stop and check grouping data and the time format argument.")
            # )
            # if(checkDate != 1) {
            #     stop('Stopped due to invalid event times.')
            # }
        }
        grouping$id <- as.character(grouping$id)
    }
    evName <- as.character(grouping$id)
    evTable <- table(evName)
    for(i in unique(evName)) {
        if(evTable[i] == 1) next
        evName[evName == i] <- paste0(i, '_',  1:evTable[i])
    }
    grouping$id <- evName
    grouping
}

# read sound acq table with minimum formatting required
readSa <- function(db) {
    con <- dbConnect(db, drv=SQLite())
    on.exit(dbDisconnect(con))
    if(!('Sound_Acquisition' %in% dbListTables(con))) {
        return(NULL)
    }
    sa <- dbReadTable(con, 'Sound_Acquisition')
    sa$Status <- str_trim(sa$Status)
    sa$SystemType <- str_trim(sa$SystemType)
    sa$UTC <- pgDateToPosix(sa$UTC)
    sa
}

nBins <- function(db) {
    evData <- getDbData(db, doSR=FALSE)
    if(nrow(evData) == 0) {
        return(0)
    }
    length(unique(evData$BinaryFile))
}

getModuleType <- function(x) {
    moduleType <- x$fileInfo$fileHeader$moduleType
    moduleType <- gsub(' ', '', moduleType)
    if(moduleType == 'SoundTrapClickDetector') {
        return('ClickDetector')
    }
    if(moduleType == 'ClickDetector' &&
       x$fileInfo$fileHeader$streamName == 'Trigger Background') {
        return('ClickDetectorTriggerBackground')
    }
    # For now, cepstrum-ing like dis
    if(moduleType == 'WhistlesMoans' &&
       grepl('cepstrum',  x$fileInfo$fileHeader$moduleName, ignore.case = TRUE)) {
        return('Cepstrum')
    }
    moduleType
}

#' @importFrom lubridate parse_date_time
#'
parseUTC <- function(x, format) {
    if(inherits(x, 'factor')) {
        x <- as.character(x)
    }
    if(is.character(x)) {
        x <- parse_date_time(x, orders=format, tz='UTC', exact=TRUE, truncated=2, quiet=TRUE)
    }
    x
}

autoMode <- function(pps, grouping) {
    if(is.data.frame(grouping) ||
       (is.character(grouping) && file.exists(grouping))) {
        return('time')
    }
    if(!is.null(getDbData(pps@db[1], doSR=FALSE)) &&
       nrow(getDbData(pps@db[1], doSR=FALSE))) {
        return('db')
    }
    if(!is.null(suppressPamWarnings(wavToGroup(pps@db[1])))) {
        return('recording')
    }
    stop('Unable to automatically detect processing mode, please ',
         'this likely means that "mode = time" is desired and no ',
         'grouping file was provided.', call.=FALSE)
}

wavToGroup <- function(db) {
    if(!file.exists(db)) {
        pamWarning('Database ', db, ' does not exist.')
        return(NULL)
    }
    con <- dbConnect(db, drv=SQLite())
    on.exit(dbDisconnect(con))
    if(!('Sound_Acquisition' %in% dbListTables(con))) {
        pamWarning('No Sound_Acquisition table in database')
        return(NULL)
    }
    sa <- dbReadTable(con, 'Sound_Acquisition')
    sa$Status <- str_trim(sa$Status)
    sa$UTC <- pgDateToPosix(sa$UTC)
    # sa <- filter(sa, .data$Status != 'Continue')
    wavCol <- findWavCol(sa)
    # Fix for merge continuos NextFile nonsense
    if(any(sa$Status == 'NextFile')) {
        whichNext <- sa$Status == 'NextFile'
        nf <- sa[whichNext, ]
        cont <- sa[sa$Status %in% c('Continue', 'Stop', 'Start'),]
        newStuff <- vector('list', length=nrow(nf))
        tempRow <- sa[FALSE,]
        tempRow[1:2, ] <- NA
        for(i in 1:nrow(nf)) {
            cBefore <- max(which(cont$Id < nf$Id[i]))
            cAfter <- min(which(cont$Id > nf$Id[i]))
            thisTemp <- tempRow
            thisTemp$Status <- c('Stop', 'Start')
            if(!is.na(wavCol)) {
                thisTemp[[wavCol]] <- c(cont[[wavCol]][cBefore],
                                        cont[[wavCol]][cAfter])
            }
            thisTemp$UTC <- cont$UTC[cAfter]
            thisTemp$Id <- nf$Id[i] + c(0, .5)
            newStuff[[i]] <- thisTemp
        }
        newStuff <- bind_rows(newStuff)
        if(is.na(wavCol)) {
            newStuff$SystemType <- unique(sa$SystemType)[1]
        }
        sa <- sa[!whichNext, ]
        sa <- rbind(sa, newStuff)
        sa <- arrange(sa, Id)
    }

    if(is.na(wavCol)) {
        pamWarning('Wav file names not saved in database, events will be labelled')
        saGrp <- select(sa, c('UTC', 'Status', 'SystemType'))
        saGrp <- filter(saGrp, .data$Status != 'Continue')
        saGrp <- distinct(arrange(saGrp, UTC))
        first <- min(which(saGrp$Status == 'Start'))
        last <- max(which(saGrp$Status == 'Stop'))
        saGrp <- saGrp[first:last,]
        alt <- saGrp$Status[1:(nrow(saGrp)-1)] != saGrp$Status[2:nrow(saGrp)]
        saGrp <- saGrp[c(TRUE, alt), ]
        saGrp$id <- rep(1:(nrow(saGrp)/2), each=2)
        saGrp$id <- paste0(gsub('\\.sqlite3', '', basename(db)), '.', saGrp$id)
        saGrp <- tidyr::spread(saGrp, 'Status', 'UTC')
        saGrp <- saGrp[!is.na(saGrp$Start) && !is.na(saGrp$Stop), ]
        if(nrow(saGrp) == 0) {
            pamWarning('Could not find appropriate start and stop times in Sound_Acquisition table')
            return(NULL)
        }

        saGrp <- data.frame(start=saGrp$Start, end=saGrp$Stop, id=saGrp$id)
    } else {
        for(i in which(is.na(sa[[wavCol]]))) {
            if(sa$Status[i] %in% c('Continue', 'Stop') &&
               i > 1 &&
               sa$Status[i-1] %in% c('Continue', 'Start') &&
               !is.na(sa[[wavCol]][i-1])) {
                sa[[wavCol]][i] <-sa[[wavCol]][i-1]
            }
            if(sa$Status[i] == 'Start' &&
               i < nrow(sa) &&
               sa$Status[i+1] %in% c('Continue', 'Stop') &&
               !is.na(sa[[wavCol]][i+1])) {
                sa[[wavCol]][i] <- sa[[wavCol]][i+1]
            }
        }
        saGrp <- select(sa, c('UTC', 'Status'))
        saGrp$id <- sa[[wavCol]]
        saGrp <- bind_rows(lapply(split(saGrp, saGrp$id), function(x) {
            thisStart <- x$UTC[grepl('Start', x$Status)]
            if(length(thisStart) == 0) {
                thisStart <- NA
            } else {
                thisStart <- max(thisStart)
            }
            thisEnd <- x$UTC[grepl('Stop', x$Status)]
            if(length(thisEnd) == 0) {
                thisEnd <- NA
            } else {
                thisEnd <- max(thisEnd)
            }
            thisId <- gsub(' ', '', unique(x$id))
            list(start = thisStart,
                 end = thisEnd,
                 id = thisId)
        }))
    }
    isNa <- is.na(saGrp$start) | is.na(saGrp$end) |
        is.infinite(saGrp$start) | is.infinite(saGrp$end)
    if(all(isNa)) {
        pamWarning('Could not find appropriate start and stop times in Sound_Acquisition table')
        return(NULL)
    }
    if(any(isNa)) {
        pamWarning('Event id(s) ', saGrp$id[isNa], ' could not get start and end times,',
                       ' they will be removed (check for missing values in Sound_Acquisition table to fix).')
        saGrp <- saGrp[!isNa,]
    }
    saGrp$db <- db
    saGrp
}

globalVariables('Id')

parseDglLabel <- function(x, colnames) {
    standardCols <- ppVars()$dglCols
    newCols <- colnames[!(colnames %in% standardCols)]
    if(length(newCols) == 0) {
        return(NA)
    }
    if(!is.null(x) &&
       x %in% newCols) {
        return(x)
    }
    if(length(newCols) == 1) {
        return(newCols)
    }
    if('Text_Annotation' %in% newCols) {
        return('Text_Annotation')
    }
    isSpec <- grepl('species', newCols, ignore.case = TRUE)
    if(sum(isSpec) == 1) {
        return(newCols[isSpec])
    }
    isLabel <- grepl('label', newCols, ignore.case=TRUE)
    if(sum(isLabel) == 1) {
        return(newCols[isLabel])
    }
    isId <- grepl('id', newCols, ignore.case=TRUE)
    if(sum(isId) == 1) {
        return(newCols[isId])
    }
    newCols[1]
}
