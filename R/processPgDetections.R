#' @title Load and Process Detections from Pamguard
#'
#' @description Loads and processes acoustic detection data that has been
#'   run through Pamguard. Uses the binary files and database(s) contained
#'   in \code{pps}, and will group your data into events by the
#'   grouping present in the 'OfflineEvents' and 'Detection Group Localiser'
#'   tables (\code{mode = 'db'}) or by the grouping specified by start and end
#'   times in the supplied \code{grouping} (\code{mode = 'time'}). Will apply
#'   all processing functions in \code{pps} to the appropriate modules
#'
#' @param pps a \linkS4class{PAMpalSettings} object containing the databases,
#'   binaries, and functions to use for processing data. See
#'   \code{\link[PAMpal]{PAMpalSettings}}. Can also be an \linkS4class{AcousticStudy}
#'   object, in which case the \code{pps} slot will be used.
#' @param mode selector for how to organize your data in to events. \code{db}
#'   will organize by events based on tables in the databases, and \code{time}
#'   will organize into events based on timestamps provided in \code{grouping}.
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
#'
#' @importFrom PamBinaries loadPamguardBinaryFile
#' @importFrom PAMmisc squishList
#' @importFrom RSQLite dbConnect dbListTables dbReadTable dbDisconnect SQLite
#' @importFrom stringr str_trim
#' @importFrom tcltk tk_choose.files
#' @importFrom purrr transpose
#' @import dplyr
#' @export
#'
processPgDetections <- function(pps, mode = c('db', 'time'), id=NULL, grouping=NULL,
                                format='%Y-%m-%d %H:%M:%OS', progress=TRUE, verbose=TRUE, ...) {
    if(identical(mode, c('db', 'time'))) {
        if(is.data.frame(grouping) ||
           (is.character(grouping) && file.exists(grouping))) {
            mode <- 'time'
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
                     'db' = processPgDetectionsDb(pps=pps, grouping=grouping, id=id,
                                                  progress=progress, ...),
                     'time' = processPgDetectionsTime(pps=pps, grouping=grouping, format=format, id=id,
                                                      progress=progress)
    )
    checkStudy(result)
    result
}

# ---- separate methods ----

#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom readr read_csv cols col_character
#'
processPgDetectionsTime <- function(pps, grouping=NULL, format='%Y-%m-%d %H:%M:%OS', id=NULL,
                                    progress=progress) {
    # start with checking grouping - parse csv if missing or provided as character and fmt times
    grouping <- checkGrouping(grouping, format)
    # this is a flag to see if any manual entries happened to grouping
    editGroup <- FALSE
    binList <- pps@binaries$list
    binFuns <- pps@functions
    allDbs <- pps@db

    # Check for what DB shit should be associated with, get full list of SA data
    # first, gonna match event times to that since its roughly the times assoicated
    # with a database
    saList <- lapply(allDbs, readSa)
    names(saList) <- allDbs

    if(!('db' %in% colnames(grouping))) {
        grouping$db <- NA_character_
    }
    # if they are there and are valid, assume they assigned
    dbToAssign <- which(!file.exists(grouping$db))
    # match db to events
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
    saByDb <- lapply(saList, function(x) unique(x$sampleRate))
    for(d in 1:nrow(grouping)) {
        if(is.na(grouping$db[d])) next
        grouping$sr[d] <- saByDb[grouping$db[d]]
    }

    # were gonna match SR by database, only need manual input if we have any
    # missing DBs
    if(any(is.na(grouping$sr))) {
        editGroup <- TRUE
        sr <- readline(prompt =
                           paste0('Not all events have a database associated with them, ',
                                  'what sample rate should be used for these events?'))
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
    modList <- c('ClickDetector', 'WhistlesMoans', 'Cepstrum')
    modWarn <- c(FALSE, FALSE, FALSE)
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
                warning('No functions for processing Module Type: ', modType, call.=FALSE)
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
            binTimes <- dplyr::bind_rows(lapply(thisBin$data, function(x) {
                list(UID = x$UID, UTC = x$date)
            }))
            binTimes$UTC <- convertPgDate(binTimes$UTC)
            binTimes <- matchSR(binTimes, thisSa)
            for(i in seq_along(thisBin$data)) {
                thisBin$data[[i]]$sr <- binTimes$sampleRate[i]
            }
        }
        thisBinData <- calculateModuleData(thisBin, binFuns)
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

    colsToDrop <- c('Id', 'comment', 'sampleRate', 'detectorName', 'parentUID',
                    'sr', 'callType', 'newUID')
    names(acousticEvents) <- evName

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
            warning('No detections in Event ', names(acousticEvents)[i], call.=FALSE)
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
            filtSa <- filter(filtSa, filtSa$UTC <= grouping$end[i], filtSa$UTC >= grouping$start[i])
            thisSource <- unique(filtSa$SystemType)
        }

        acousticEvents[[i]] <-
            AcousticEvent(id=evName[i], detectors = thisData, settings = list(sr = thisSr, source = thisSource),
                          files = list(binaries=binariesUsed, db=grouping$db[i], calibration=calibrationUsed))
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
                           files = list(db=allDbs, binaries=allBins),
                           ancillary = list(grouping=grouping))
    on.exit() # this cancels the on.exit 'save my grouping' call that is there if you crash
    study
}

#'
processPgDetectionsDb <- function(pps, grouping=c('event', 'detGroup'), id=NULL,
                                  progress=TRUE, ...) {
    allDb <- pps@db
    # awk diff init values between modes have to reset this here
    if(is.null(grouping)) {
        grouping <- c('event', 'detGroup')
    }

    nBin <- sum(sapply(allDb, nBins))
    if(progress) {
        cat('Processing databases... \n')
        pb <- txtProgressBar(min=0, max=nBin, style=3)
    }
    binNo <- 1
    modList <- c('ClickDetector', 'WhistlesMoans', 'Cepstrum')
    modWarn <- c(FALSE, FALSE, FALSE)
    names(modWarn) <- modList
    allAcEv <- lapply(allDb, function(db) {
        tryCatch({
            missBins <- character(0)
            failBin <- 'No file processed'
            binList <- pps@binaries$list
            binFuns <- pps@functions
            dbData <- getDbData(db, grouping, ...)
            if(is.null(dbData) ||
               nrow(dbData) == 0) {
                warning('No detections found in database ',
                        basename(db), '.', call.=FALSE)
                # setTxtProgressBar(pb, value = evNo)
                # evNo <- evNo + 1
                return(NULL)
            }
            thisSr <- unique(dbData$sampleRate)
            if(length(thisSr) > 1) {
                warning('More than 1 sample rate found in database ',
                        basename(db),'.', call.=FALSE)
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
                            warning('No functions for processing Module Type: ', modType, call.=FALSE)
                        }
                    }
                    binData <- calculateModuleData(thisBin, binFuns)
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
                warning('Could not find the matching binary files for binaries ',
                        printN(missBins, 3),
                        ' in database ', basename(db), call.=FALSE)
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
                lapply(dbData, function(x) split(x, x$parentUID)),
                .names = unique(unlist(sapply(dbData, function(x) x$parentUID)))
            )

            # Should this function store the event ID? Right now its just the name
            # in the list, but is this reliable? Probably not
            colsToDrop <- c('Id', 'comment', 'sampleRate', 'detectorName', 'parentUID',
                            'sr', 'callType', 'newUID')
            acousticEvents <- lapply(dbData, function(ev) {
                ev <- ev[sapply(ev, function(x) !is.null(x))]
                binariesUsed <- sapply(ev, function(x) unique(x$BinaryFile)) %>%
                    unlist(recursive = FALSE) %>% unique()
                binariesUsed <- unlist(sapply(binariesUsed, function(x) grep(x, binList, value=TRUE, fixed=TRUE), USE.NAMES = FALSE))
                evId <- paste0(gsub('\\.sqlite3', '', basename(db)), '.', unique(ev[[1]]$parentUID))
                ev <- lapply(ev, function(x) {
                    x$BinaryFile <- basename(x$BinaryFile)
                    thisType <- unique(x$callType)
                    x <- dropCols(x, colsToDrop)
                    attr(x, 'calltype') <- thisType
                    x
                })
                AcousticEvent(id = evId, detectors = ev, settings = list(sr = thisSr, source=thisSource),
                              files = list(binaries=binariesUsed, db=db, calibration=calibrationUsed))
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
    AcousticStudy(id=id, events = allAcEv, pps = pps,
                  files = list(db=allDbs, binaries=allBins))
}

# ---- not exported helpers ----
getDbData <- function(db, grouping=c('event', 'detGroup'), label=NULL) {
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
            suppressWarnings(
                bind_rows(
                    lapply(grouping, function(x) {
                        getDbData(db, x)
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
               eventColumns <- c('UID', label)
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
               if(is.null(label)) {
                   label <- 'Text_Annotation'
               }
               eventColumns <- c('UID', label)
               evName <- 'DGL'
           },
           {
               stop("I don't know how to group by ", grouping, '.\n', call.=FALSE)
           }
    )

    if(length(detTables)==0 ||
       length(eventTables)==0) {
        warning('Could not find event tables for grouping method "', grouping,
                '" in database ', basename(db), call.=FALSE)
        return(NULL)
    }
    allDetections <- bind_rows(
        lapply(detTables, function(table) {
            dbReadTable(con, table)
        })
    )
    if(nrow(allDetections)==0) {
        warning('No detections found for grouping method "', grouping,
                '" in database ', basename(db), call.=FALSE)
        return(NULL)
    }

    allEvents <- bind_rows(
        lapply(eventTables, function(table) {
            dbReadTable(con, table)
        })
    )
    if(nrow(allEvents)==0) {
        warning('No events found for grouping method "', grouping,
                '" in database ', basename(db), call.=FALSE)
        return(NULL)
    }

    eventColumns <- eventColumns[eventColumns %in% colnames(allEvents)]
    allEvents <- select_(allEvents, .dots=eventColumns)

    # Do i want all detections in clicks, or only all in events?
    # left_join all det, inner_join ev only
    if(!('UID' %in% names(allEvents)) ||
       !('parentUID' %in% names(allDetections))) {
        message('UID and parentUID columns not found in database ', basename(db),
                ', these are required to process data. Please upgrade to Pamguard 2.0+.')
        return(NULL)
    }

    allDetections <- inner_join(
        allDetections, allEvents, by=c('parentUID'='UID')
    )
    if(!('newUID' %in% colnames(allDetections))) {
        allDetections$newUID <- -1
    }
    allDetections <- allDetections %>%
        mutate(BinaryFile = str_trim(.data$BinaryFile),
               # UTC = as.POSIXct(as.character(UTC), format='%Y-%m-%d %H:%M:%OS', tz='UTC')) %>%
               UTC = pgDateToPosix(.data$UTC)) %>%
        select_(.dots=unique(c(eventColumns, 'UTC', 'UID', 'parentUID', 'BinaryFile', 'newUID')))

    # rename column to use as label - standardize across event group types
    colnames(allDetections)[which(colnames(allDetections)==label)] <- 'eventLabel'

    allDetections <- matchSR(allDetections, db, extraCols=c('SystemType'))

    # apply str_trim to all character columns
    whichChar <- which(sapply(allDetections, function(x) 'character' %in% class(x)))
    for(i in whichChar) {
        allDetections[, i] <- str_trim(allDetections[, i])
    }
    allDetections <- select(allDetections, -.data$UTC)
    allDetections$UID <- as.character(allDetections$UID)
    allDetections$newUID <- as.character(allDetections$newUID)
    allDetections$parentUID <- paste0(evName, allDetections$parentUID)
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
        warning(paste0('UID(s) ', printN(dbData$UID[!dbData$matched], n=6),
                       ' are in databases ', dbName, ' but not in binary file ', binFile), call.=FALSE)
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
        grouping <- read_csv(grouping, col_types = cols(.default=col_character()))
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
    sa <- dbReadTable(con, 'Sound_Acquisition')
    sa$Status <- str_trim(sa$Status)
    sa$SystemType <- str_trim(sa$SystemType)
    sa$UTC <- pgDateToPosix(sa$UTC)
    sa
}

nBins <- function(db) {
    evData <- getDbData(db)
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
