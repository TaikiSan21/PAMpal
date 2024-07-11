#' @title Match time-based data to PAMpal objects
#'
#' @description Match any time-based data (dataframe with a \code{UTC} column)
#'   to an AcousticStudy or AcousticEvent object
#'
#' @param x \linkS4class{AcousticStudy} or \linkS4class{AcousticEvent} object
#'   to match data to, or a dataframe for \code{timeJoin}
#' @param data a data frame to match to data to \code{x}.
#'   Must have column \code{UTC}, and optionally a column \code{db} if subsets
#'   of data should be matched only to parts of \code{x} with that database. All
#'   other columns will be considered variables to add to \code{x}
#' @param mode one of "event" or "detection". "event" will match one set of
#'   variables per event, stored in the "measures" for that event. "detection"
#'   will match variables to every detection.
#' @param thresh maximum time apart in seconds for matching to
#'   data, if the closest value is more than \code{thresh} apart then the
#'   variable values will be set to \code{NA}
#' @param interpolate logical flag whether or not to interpolate between points
#'   in \code{data} or just matched to nearest time
#' @param replace one of \code{TRUE}, \code{FALSE}, or \code{NA}. If \code{TRUE},
#'   all existing values with the same name as columns in \code{data} will be
#'   replaced. If \code{FALSE} no replacement occurs. If \code{NA} only values
#'   which are \code{NA} will be replaced with new values
#' @param keepDiff logical flag to keep time difference data
#'
#' @details This function lets you match any arbitrary data to a PAMpal object
#'   as long as it has a time associated with it. Data will be attached to
#'   detector dataframes for \code{mode='detection'} or to the event "measures"
#'   location for \code{mode='event'} (this is where \link{calculateICI} and
#'   \link{matchEnvData,AcousticStudy-method} store their event data). These can be accessed with the
#'   \link{getMeasures} function and are also exported in the various "getXXX"
#'   functions (\link{getClickData} etc.) if \code{measures=TRUE} (default).
#'
#'   All columns in the provided \code{data} object will be treated as variables
#'   to add, with a few exceptions. There are a few reserved column names used by
#'   PAMpal that cannot be overridden (e.g. UID, eventId, species). Also any columns
#'   already existing in the PAMpal data will not be overridden unless \code{replace}
#'   is not \code{FALSE}. The column names in \code{data} will be used as the names
#'   for the added variables, so care should be taken to ensure these are informative
#'   enough for future use.
#'
#' @return the same data as \code{x}, with data added from \code{data}
#'
#' @author Taiki Sakai \email{taiki.sakai@@noaa.gov}
#'
#' @examples
#'
#' data(exStudy)
#' addData <- data.frame(UTC = as.POSIXct('2018-03-20 15:25:10', tz='UTC'),
#'                       newVariable = 26)
#' data <- matchTimeData(exStudy, addData, mode='detection')
#' getClickData(data)
#' data <- matchTimeData(exStudy, addData, mode='event')
#' getMeasures(data)
#'
#' @export
#'
matchTimeData <- function(x, data, mode=c('event', 'detection'), thresh=Inf, interpolate=TRUE, replace=FALSE, keepDiff=FALSE) {
    if(is.AcousticStudy(x)) {
        dbFilt <- FALSE
        if('db' %in% colnames(data)) {
            xDb <- files(x)$db
            if(!any(basename(xDb) %in% basename(unique(data$db)))) {
                stop('"db" column in "data" does not match any databases in "x".',
                     ' Check spelling or remove "db" column to continue.')
            }
            data <- split(data, basename(data$db))
            # make fake missing data so that we can fill with NA so all events
            # keep same dimensions. other wise getXXX functions will break
            missing <- basename(xDb)[!basename(xDb) %in% names(data)]
            if(length(missing) > 0) {
                fakeData <- data[[1]][1, ]
                fakeData[1, ] <- NA
                fakeData$UTC <- data[[1]]$UTC[1]
                for(m in missing) {
                    fakeData$db <- m
                    data[[m]] <- fakeData
                }
            }
            dbFilt <- TRUE
        }
        for(e in seq_along(events(x))) {
            if(dbFilt) {
                x[[e]] <- matchTimeData(x[[e]], data[[basename(files(x[[e]])$db)]], mode, thresh, interpolate, replace)
            } else {
                x[[e]] <- matchTimeData(x[[e]], data, mode, thresh, interpolate, replace)
            }
        }
        x <- .addPamWarning(x)
        return(x)
    }
    if(!is.AcousticEvent(x)) {
        stop('"x" must be AcousticStudy or AcousticEvent')
    }
    switch(
        match.arg(mode),
        'event' = {
            eventData <- getEventStart(x)
            eventData <- timeJoin(eventData, data, thresh=thresh, interpolate=interpolate, replace=replace)
            if(eventData$timeDiff > thresh) {
                pamWarning('Event ', id(x), ' exceeded time threshold, setting values',
                           ' to NA. (Average ', round(mean(eventData$timeDiff/3600, na.rm=TRUE),1), ' hours apart)',
                           which=sys.nframe()-1)
            }
            toDrop <- c('UTC', 'eventId')
            if(isFALSE(keepDiff)) {
                toDrop <- c(toDrop, 'timeDiff')
            }
            eventData <- as.list(dropCols(eventData, toDrop))
            # REPLACE DOES NOTHING HERE
            oldMeas <- ancillary(x)$measures
            if(!is.null(oldMeas)) {
                eventData <- safeListAdd(oldMeas, eventData, replace=replace)
            }
            ancillary(x)$measures <- eventData
        },
        'detection' = {
            diffs <- numeric(0)
            for(d in seq_along(detectors(x))) {
                calltype <- attr(x[[d]], 'calltype')
                thisDet <- x[[d]]
                oldCols <- colnames(thisDet)
                thisDet <- timeJoin(thisDet, data, thresh=thresh, interpolate=interpolate, replace=replace)
                newNa <- sapply(colnames(thisDet)[!colnames(thisDet) %in% oldCols], function(c) {
                    anyNA(thisDet[[c]])
                })
                newNa <- map(colnames(thisDet)[!colnames(thisDet) %in% oldCols], function(c) {
                             is.na(thisDet[[c]])
                }) %>%
                    reduce(`|`)
                if(any(newNa)) {
                    diffs <- c(diffs, thisDet$timeDiff[newNa])
                }
                if(isFALSE(keepDiff)) {
                    thisDet$timeDiff <- NULL
                }
                attr(thisDet, 'calltype') <- calltype
                x[[d]] <- thisDet
                # x[[d]] <- timeJoin(x[[d]], data, thresh=thresh, interpolate=interpolate, replace=replace)

                # attr(x[[d]], 'calltype') <- calltype
            }
            if(length(diffs) > 0) {
                pamWarning(length(diffs), ' matches in event ', id(x), ' exceeded time threshold, setting',
                           ' to NA. (Average ', round(mean(diffs/3600, na.rm=TRUE),1), ' hours apart)',
                           which=sys.nframe()-1)
            }
        }
    )
    x
}

#' @rdname matchTimeData
#'
#' @param y dataframe to join to \code{x}
#'
#' @export
#'
timeJoin <- function(x, y, thresh=Inf, interpolate=TRUE, replace=FALSE) {
    #add reserved column names like eventId, db
    if(is.null(x) ||
       nrow(x) == 0) {
        return(x)
    }
    x[['timeDiff']] <- NULL
    reservedCols <- c('eventId', 'db', 'UID', 'eventLabel', 'BinaryFile', 'species', 'Channel')
    if(isFALSE(replace)) {
        addCols <- colnames(y)[!colnames(y) %in% colnames(x)]
        oldCols <- colnames(x)
    } else if(isTRUE(replace)) {
        addCols <- colnames(y)[colnames(y) != 'UTC']
        oldCols <- colnames(x)[!colnames(x) %in% colnames(y)]
        oldCols <- c('UTC', oldCols)
    } else if(is.na(replace)) {
        overlaps <- (colnames(y) %in% colnames(x)) &
            colnames(y) != 'UTC'
        if(any(overlaps)) {
            overlaps <- colnames(y)[overlaps]
            hasNa <- sapply(overlaps, function(c) {
                anyNA(x[[c]])
            })
            hasNa <- overlaps[hasNa]
            for(n in hasNa) {
                isNa <- is.na(x[[n]])
                newMatch <- timeJoin(x[c('UTC', n)], y[c('UTC', n)], thresh=thresh, interpolate=interpolate, replace=TRUE)
                x[[n]][isNa] <- newMatch[[n]][isNa]
            }
        }
        addCols <- colnames(y)[!colnames(y) %in% colnames(x)]
        oldCols <- colnames(x)
    }
    addCols <- addCols[!addCols %in% reservedCols]
    if(length(addCols) == 0) {
        return(x)
    }
    y <- y[c('UTC', addCols)]
    x <- x[oldCols]
    x$xTime <- x$UTC
    y$yTime <- y$UTC
    setDT(x)
    setDT(y)
    setkeyv(x, 'UTC')
    setkeyv(y, 'UTC')
    result <- y[x, roll='nearest']
    result$timeDiff <- abs(as.numeric(result$xTime) - as.numeric(result$yTime))
    tooFar <- result$timeDiff > thresh
    result[tooFar, (addCols) := NA]
    if(nrow(y) == 1) {
        interpolate <- FALSE
    }
    if(interpolate) {
        for(c in addCols) {
            interpVal <- approx(x=y$UTC, y=y[[c]], xout=result$xTime[!tooFar], rule=2)$y
            result[[c]][!tooFar] <- interpVal
        }
    }
    result[, c('xTime', 'yTime') := NULL]
    setDF(x)
    setDF(y)
    setDF(result)
    result[c(oldCols, addCols, 'timeDiff')]
}
