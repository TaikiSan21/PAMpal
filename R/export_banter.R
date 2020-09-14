#' @title Export a List of AcousticEvent Ojbects for Banter Model
#'
#' @description Formats a list of AcousticEvent objects into the structure
#'   needed to run a banter model.
#'
#' @param x a \linkS4class{AcousticStudy} object or
#'   a list of \linkS4class{AcousticEvent} objects
#' @param dropVars a vector of the names of any variables to remove
#' @param dropSpecies a vector of the names of any species to exclude
#' @param training logical flag whether or not this will be used as a
#'   training data set, or a value between 0 and 1 specifying what percent
#'   of the data should be used for training (with the rest set aside
#'   for testing). If TRUE or greater than 0, must contain species ID.
#'   NOTE: if value is not 0, 1, \code{TRUE}, or \code{FALSE}, output
#'   will be further split into \code{training} and \code{test} items
#'   within the list output
#'
#' @return a list with three items, \code{events}, \code{detectors}, and
#'   \code{na}. If value of \code{training} is not 0, 1, \code{TRUE}, or
#'   \code{FALSE}, output will be split into \code{training} and
#'   \code{test} lists that contain \code{events} and \code{detectors}.
#'   \code{events} is a dataframe with two columns. \code{event.id} is a
#'   unique identifier for each event, taken from the names of the event
#'   list. \code{species} is the species classification, taken from the
#'   \code{species} slot labelled \code{id}. \code{detectors} is a list
#'   of data frames containing all the detections and measurements. There is
#'   one list for each unique detector type found in the \code{detectors} slots
#'   of \code{x}. The data frames will only have columns with class
#'   \code{numeric}, \code{integer}, \code{factor}, or \code{logical}, and
#'   will also have columns named \code{UID}, \code{Id}, \code{parentUID},
#'   \code{sampleRate}, \code{Channel}, \code{angle}, and \code{angleError},
#'   removed so that these are not treated as parameters for the banter random
#'   forest model. The dataframes will also have columns \code{event.id} and
#'   \code{call.id} added. \code{na} contains the UIDs and Binary File names
#'   for any detections that had NA values. These cannot be used in the
#'   random forest model and are removed from the exported dataset.
#'
#' @author Taiki Sakai \email{taiki.sakai@@noaa.gov}
#'
#' @importFrom PAMmisc squishList
#' @importFrom dplyr distinct n group_by n_distinct summarise rename left_join bind_rows filter
#' @importFrom knitr kable
#' @importFrom purrr map reduce transpose
#' @export
#'
export_banter <- function(x, dropVars=NULL, dropSpecies=NULL, training=TRUE) {
    if(is.list(x)) {
        whichStudy <- sapply(x, is.AcousticStudy) # two cases if its a list
        whichEvent <- sapply(x, is.AcousticEvent)
        # NEED TO UNLIST IF ITS A STUDY BUT THEN PROB SHOULDNT DO ALL IN ONE THING WHATEVER BRO
        x <- c(x[whichEvent],
               unlist(sapply(x[whichStudy], events)))
        # x <- sapply(x, function(e) {
        #     if(is.AcousticStudy(e)) {
        #         return(events(e))
        #     } else if(is.AcousticEvent(e)) {
        #         return(e)
        #     } else {
        #         warning('Some inputs were not an AcousticEvent or AcousticStudy.')
        #         return(NULL)
        #     }
        # })
    }
    if(is.AcousticStudy(x)) {
        x <- events(x)
    }
    if(is.AcousticEvent(x)) {
        x <- list(x)
        names(x) <- x[[1]]@id
    }
    # From here we need x to be list of AcEv
    numDets <- sapply(x, function(e) length(detectors(e)))
    noDets <- (numDets == 0)
    if(all(noDets)) {
        stop('No events had any detections.')
    }
    if(any(noDets)) {
        warning(sum(noDets), ' events have 0 detections, they will be removed.\n', call. = FALSE)
        x <- x[!noDets]
    }
    sp <- sapply(x, function(e) species(e)$id)
    sp[is.null(sp)] <- NA_character_
    spNa <- sapply(sp, is.na)
    if(training && any(spNa)) {
        warning('Events ', paste(names(x)[which(spNa)], collapse=', '),
             ' do not have a species ID. Data can only be used for prediction, not model training.', call. = FALSE)
        training <- FALSE
    }

    detNA <- data.frame(UID = character(0), BinaryFile = character(0),
                        event = character(0), detector = character(0), stringsAsFactors = FALSE)
    evName <- sapply(x, id)
    if(!(length(unique(evName)) == length(evName))) {
        warning('Duplicate event names found, these must be unique for BANTER. Adding numbers to event names.', call. = FALSE)
        for(i in unique(evName)) {
            evName[evName == i] <- paste0(i, 1:(sum(evName == i)))
        }

    }
    names(x) <- evName
    events <- data.frame(event.id = names(x),
                         stringsAsFactors = FALSE)

    # check if enough events to make train data, remove any we dont have for
    if(training) {
        events$species <- sp
        nSpecies <- table(events$species)
        # need min 2 to train banter, 3 if also making test
        if(training == 1) {
            noTrainSpecies <- names(nSpecies)[nSpecies < 2]
            noTrainSpecies <- noTrainSpecies[!(noTrainSpecies %in% dropSpecies)]
            if(length(noTrainSpecies) > 0) {
                warning('Species ', paste0(noTrainSpecies, collapse=', '),
                        ' do not have enough events to train a banter model (min 2),',
                        ' these will be removed.', call. = FALSE)
                dropSpecies <- c(dropSpecies, noTrainSpecies)
            }
        } else {
            noTrainSpecies <- names(nSpecies)[nSpecies < 3]
            noTrainSpecies <- noTrainSpecies[!(noTrainSpecies %in% dropSpecies)]
            if(length(noTrainSpecies) > 0) {
                warning('Species ', paste0(noTrainSpecies, collapse=', '),
                        ' do not have enough events to train a banter model (min 2',
                        ' to train plus 1 for testing), these will be removed.', call. = FALSE)
                dropSpecies <- c(dropSpecies, noTrainSpecies)
            }
        }
    }

    if(!is.null(dropSpecies)) {
        toDrop <- sp %in% dropSpecies
        sp <- sp[!toDrop]
        events <- events[!toDrop, , drop=FALSE] # dont coerce to vector ever
        x <- x[!toDrop]
    }
    if(nrow(events) == 0) {
        stop('No events left with any detections after removing species.')
    }
    # check if any event level measures are present in all data
    # or should i get all and fill NA, let banter deal with the NAs and warn you?
    # if any present in all get them and cbind that ish to your events
    allMeasures <- reduce(
        lapply(x, function(e) names(ancillary(e)$measures)),
        intersect)
    # if(length(measureNames) == 0) {
    #     allMeasures <- NULL
    # } else {
    #     allMeasures <- reduce(measureNames, intersect)
    # }
    if(length(allMeasures) > 0 ) {

        measureData <- bind_rows(lapply(x[events[['event.id']]], function(e) {
            ancillary(e)$measures[allMeasures]
        }))
        events <- cbind(events, measureData)
        cat('Found ', length(allMeasures), ' event level measures that were present',
            ' in all events, adding these to your event data.\n', sep='')
    }

    for(e in seq_along(x)) {
        thisEv <- x[[e]]
        for(d in seq_along(detectors(thisEv))) {
            thisDet <- detectors(thisEv)[[d]]
            if(is.null(thisDet)) next

            thisDet$event.id <- names(x)[e]
            thisDet$call.id <- paste0(names(x)[e], thisDet$UID)
            if('Channel' %in% colnames(thisDet)) {
                thisDet$call.id <- paste0('C', thisDet$Channel, thisDet$call.id)
            }
            colsToDrop <- c('UID', 'Id', 'parentUID', 'sampleRate', 'Channel',
                            'angle', 'angleError', 'peakTime')
            colsToDrop <- unique(c(colsToDrop, dropVars))
            useCols <- lapply(thisDet, class) %in% c('numeric', 'integer', 'factor', 'logical') &
                !(colnames(thisDet) %in% colsToDrop) |
                colnames(thisDet) %in% c('event.id', 'call.id')

            whereNA <- sapply(thisDet[, useCols], is.na)
            if(nrow(thisDet) == 1) {
                whereNA <- matrix(whereNA, nrow=1)
            }
            naRow <- apply(whereNA, 1, any)
            if(length(naRow) != nrow(thisDet)) browser()
            thisNA <- thisDet[naRow, c('UID', 'BinaryFile')]
            if(nrow(thisNA) > 0) {
                thisNA$event <- names(x)[e]
                thisNA$detector <- names(detectors(thisEv))[d]
                thisNA$measureNames <- apply(whereNA, 1, function(x) paste0(colnames(thisDet[, useCols])[x], collapse=', '))[naRow]
            }
            detNA <- rbind(detNA, thisNA)
            detectors(thisEv)[[d]] <- thisDet[!naRow, useCols]
        }
        x[[e]] <- thisEv
    }

    dets <- lapply(x, function(e) {
        tmpDet <- detectors(e)
        if(length(tmpDet) == 0) return(NULL)
        tmpDet[sapply(tmpDet, function(y) !is.null(y) && ncol(y) > 2)]
    })
    names(dets) <- NULL
    dets <- squishList(unlist(dets, recursive = FALSE))
    dets <- lapply(dets, distinct)
    if(nrow(detNA) > 0) {
        warning('Removing ', nrow(detNA), ' NA values, to see affected UID(s) and ',
        'BinaryFile(s) check the "na" item in list output.', call. = FALSE)
    }
    # From here train/test split and report
    # result <- list(events=events, detectors=dets, na=detNA)
    result <- bntSplit(list(events=events, detectors=dets), training)
    # first two cases return $events and $detectors
    if(training == 0) {
        cat('\nCreated data for ', nrow(result$events), ' events with ',
        sum(sapply(result$detectors, nrow)), ' total detections.', sep='')
    } else if(training==1) {
        cat('\nCreated data for ', nrow(result$events), ' events with ',
            sum(sapply(result$detectors, nrow)), ' total detections', sep='')
        cat(' and ', length(unique(sp)),
            ' unique species: ', paste0(unique(sp), collapse =', '), '.',
            '\nRe-run with dropSpecies argument if any of these are not desired.', sep='')
        print(kable(bntSummaryTable(result)))
    } else if(training > 0) { # train-test split case has $train$events and $test$events
        cat('\nCreated training data for ', nrow(result$train$events), ' events with ',
            sum(sapply(result$train$detectors, nrow)), ' total detections', sep='')
        cat(' and ', length(unique(sp)),
            ' unique species: ', paste0(unique(sp), collapse =', '), '.', sep='')

        print(kable(bntSummaryTable(result$train)))
        cat('\nCreated test data for ', nrow(result$test$events), ' events with ',
            sum(sapply(result$test$detectors, nrow)), ' total detections', sep='')
        print(kable(bntSummaryTable(result$test)))
        cat('\nRe-run with dropSpecies argument if any of these are not desired.', sep='')
    }
    cat('\n')
    # list(events=events, detectors=dets, na=detNA)
    result$na <- detNA
    result
}

# make a summary table from export_banter output
bntSummaryTable <- function(x) {
    bind_rows(lapply(x$detectors, function(y) {
        tmp <- left_join(y, x$events, by='event.id')
        select(tmp, .data$event.id, .data$species)
    }), .id = 'detector') %>% group_by(species) %>%
        summarise(Events = n_distinct(.data$event.id),
                  Detections = n()) %>%
        rename(Species = species)
}

# x is result from initial export_banter with $events and $detectors for all
bntSplit <- function(x, trainFrac) {
    testFrac <- 1-trainFrac
    if(testFrac == 0 ||
       trainFrac == 0) {
        return(x)
    }
    split <- split(x$events, x$events$species)
    ev <- transpose(lapply(split, function(s) {
        if(nrow(s) < 3) {
            warning('Species ', s$species[1],
                    ' does not have enough events to train a model (min 3).', call. = FALSE)
            return(NULL)
        }
        if(nrow(s) * testFrac < 1) {
            warning('Species ', s$species[1],
                    ' does not have enough events for the chosen test proportion,',
                    ' will round up to 1 event for test set.', call. = FALSE)
        }
        testNum <- ceiling(nrow(s) * testFrac)
        testIx <- sample(1:nrow(s), testNum, replace=FALSE)
        trainIx <- (1:nrow(s))[-testIx]
        list(train=s[trainIx, ], test=s[testIx, ])
    }))
    result <- list(train = list(events = bind_rows(ev$train)),
                   test = list(events = bind_rows(ev$test)))
    # result$train$detectors <- x$detectors$event.id %in% train$event.id #psuedodfacode
    result$train$detectors <- lapply(x$detectors, function(d) {
        tmp <- filter(d, .data$event.id %in% result$train$events$event.id)
        # if(nrow(tmp) == 0) return(NULL)
        tmp
    })
    result$test$detectors <- lapply(x$detectors, function(d) {
        tmp <- filter(d, .data$event.id %in% result$test$events$event.id)
        # if(nrow(tmp) == 0) return(NULL)
        tmp
    })
    result
}
