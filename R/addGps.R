#' @title Add GPS Locations to an AcousticStudy
#'
#' @description Add GPS Lat / Long to an AcousticStudy or AcousticEvent.
#'   If GPS data is not present in any of the databases, user will
#'   interactively be asked to provide GPS data to add
#'
#' @param x data to add GPS coordinates to. Must have a column \code{UTC}, and
#'   can also have an optional column \code{Channel}
#' @param gps a data frame of GPS coordinates to match to data from \code{x}.
#'   Must have columns \code{UTC}, \code{Latitude}, \code{Longitude}, and
#'   optionally \code{Channel}. If not provided and \code{x} is an
#'   \linkS4class{AcousticEvent} or \linkS4class{AcousticStudy} object, then
#'   the gps data will be read from the databases contained in the \code{files}
#'   slot of \code{x}
#' @param thresh maximum time apart in seconds for matching GPS coordinates to
#'   data, if the closest coordinate is more than \code{thresh} apart then the
#'   Latitude and Longitude values will be set to \code{NA}
#' @param keepDiff logical flag to keep time difference column (between GPS time
#'   and data time)
#' @param \dots additional arguments for other methods
#'
#' @details Latitude and Longitude coordinates will be matched to the data
#'   by using data.tables rolling join with \code{roll='nearest'}. After the
#'   join is done, the time difference between the matched rows is checked
#'   and any that are greater than the set threshold are set to NA. This is
#'   done to prevent accidentally matching weird things if an incomplete set
#'   of GPS data is provided.
#'
#'   If \code{x} is an \linkS4class{AcousticEvent} or \linkS4class{AcousticStudy},
#'   then \code{gps} can be omitted and will be read from the databases contained
#'   in the \code{files} slot of \code{x}. If \code{x} is an \linkS4class{AcousticStudy},
#'   then the gps data will also be saved to the \code{gps} slot of the object, and
#'   an additional argument \code{bounds} can be provided. This is a length two vector
#'   of \code{POSIXct} class times that will bound the times of gps data to store, gps
#'   data outside this range will not be stored (to reduce the potentially very large
#'   amount of data stored in the \code{gps} slot)
#'
#' @return the same data as \code{x}, with Lat/Long data added. AcousticStudy objects
#'   will have all GPS data used added to the "gps" slot, and all AcousticEvents will
#'   have Latitude and Longitude added to all detector dataframes
#'
#' @author Taiki Sakai \email{taiki.sakai@@noaa.gov}
#'
#' @examples
#'
#' data(exStudy)
#' # need to update database file to local directory
#' db <- system.file('extdata', 'Example.sqlite3', package='PAMpal')
#' exStudy <- updateFiles(exStudy, db=db, bin=NA, verbose=FALSE)
#' exStudy <- addGps(exStudy)
#' head(gps(exStudy))
#'
#' @name addGps
#' @export
#'
setGeneric('addGps', function(x, gps=NULL, thresh = 3600, ...) standardGeneric('addGps'))

#' @rdname addGps
#' @importFrom data.table data.table setkeyv key setDT setDF
#' @importFrom dplyr mutate select
#' @importFrom utils globalVariables
#' @importFrom stats approxfun
#' @importFrom geosphere distGeo
#' @export
#'
setMethod('addGps', 'data.frame', function(x, gps, thresh = 3600, keepDiff=FALSE, ...) {
    # Check for right columns and proper types
    needCols <- c('UTC', 'Latitude', 'Longitude')
    missingCols <- needCols[!(needCols %in% colnames(gps))]
    if(length(missingCols) > 0) {
        pamWarning('Gps data needs column(s) named ', paste(missingCols, collapse = ', '))
        return(x)
    }
    if(!('UTC' %in% colnames(x))) {
        pamWarning('Data needs column UTC.')
        return(x)
    }
    if(!('POSIXct' %in% class(x$UTC))) x$UTC <- pgDateToPosix(x$UTC)
    if(!('POSIXct' %in% class(gps$UTC))) gps$UTC <- pgDateToPosix(gps$UTC)
    # code for just interping:
    # can just make functions once at start and pass along
    gps <- distinct(gps)
    gpsRange <- range(gps$UTC, na.rm=TRUE)
    interpMethod <- ifelse(nrow(gps) > 1, 'linear', 'constant')
    latFun <- approxfun(x=gps$UTC, y=gps$Latitude, rule=2, method=interpMethod, ties='ordered')
    lonFun <- approxfun(x=gps$UTC, y=gps$Longitude, rule=2, method=interpMethod, ties='ordered')
    # then use these for interpd results
    # x$Latitude <- latFun(x$UTC)
    # x$Longitude <- lonFun(x$UTC)
    # problem with this is that we get nothing about really far apart matches
    # dummies for calculating time difference for threshold check later
    thisType <- attr(x, 'calltype')
    # x$dataTime <- x$UTC
    # gps$gpsTime <- gps$UTC
    x <- dropCols(x, c('Latitude', 'Longitude'))
    setDT(x)
    # x[, dataTime := UTC]
    x$dataTime <- x$UTC
    gps <- checkGpsKey(gps)
    
    if('Channel' %in% colnames(gps)) {
        gpsCols <- c('Channel', needCols)
        # gps <- gps[, c('Channel', 'gpsTime', needCols), with=FALSE]
    } else {
        gpsCols <- needCols
        # gps <- gps[, c('gpsTime', needCols), with=FALSE]
    }
    
    gps <- gps[, gpsCols, with=FALSE]
    # gps[, gpsTime := UTC]
    gps$gpsTime <- gps$UTC
    # set keys for the join, using 'Channel' if its present
    if('Channel' %in% colnames(x) &&
       'Channel' %in% colnames(gps)) {
        setkeyv(x, c('Channel', 'UTC'))
    } else {
        setkeyv(x, 'UTC')
        setkeyv(gps, 'UTC') # removing channel key from gps if its there i guess
    }
    result <- gps[x, roll='nearest']
    result$closeLat <- result$Latitude
    result$closeLon <- result$Longitude
    # result[, closeLat := Latitude]
    # result[, closeLon := Longitude]
    tooFar <- abs(as.numeric(difftime(result$dataTime, result$gpsTime, units='secs'))) > thresh
    if(keepDiff) {
        result$timeDiff <- abs(as.numeric(difftime(result$dataTime, result$gpsTime, units='secs')))
    }
    result[tooFar, c('Latitude', 'Longitude') := NA]
    # result[abs(as.numeric(difftime(dataTime, gpsTime, units='secs'))) > thresh, c('Latitude', 'Longitude') := NA]
    # result$Latitude[tooFar] <- NA
    # result$Longitude[tooFar] <- NA
    result$Latitude[!tooFar] <- latFun(result$dataTime[!tooFar])
    result$Longitude[!tooFar] <- lonFun(result$dataTime[!tooFar])
    result$gpsUncertainty <- NA
    result$gpsUncertainty[!tooFar] <- 
        distGeo(matrix(c(result$Longitude[!tooFar], result$Latitude[!tooFar]), ncol=2),
                           matrix(c(result$closeLon[!tooFar], result$closeLat[!tooFar]), ncol=2))
    result[, UTC := dataTime]
    # if were interp'ing out of bounds this ends up at 0 which is incorrect
    if(interpMethod == 'linear') {
        result$gpsUncertainty[result$UTC <= gpsRange[1]] <- NA
        result$gpsUncertainty[result$UTC >= gpsRange[2]] <- NA
    }
    # result$UTC <- result$dataTime
    if(any(is.na(result$Longitude))) {
        avgTime <- mean(abs(as.numeric(difftime(result$dataTime[tooFar], result$gpsTime[tooFar], units='hours'))))
        warning('Some GPS coordinate matches exceeded time threshold, setting',
                ' value to NA. (Average ', round(avgTime,1), ' hours apart)')
    }
    result[, c('gpsTime', 'dataTime', 'closeLat', 'closeLon') := NULL]
    attr(result, 'calltype') <- thisType
    setDF(result)
    result
})

#' @rdname addGps
#' @export
#'
setMethod('addGps', signature(x='AcousticEvent'), function(x, gps=NULL, thresh = 3600, ...) {
    if(is.null(gps)) {
        gps <- rbindlist(lapply(files(x)$db, function(db) {
            # gpsFromDb(db, extraCols=c('Speed', 'Heading', 'MagneticVariation', 'db'), ...)
            gpsFromDb(db, extraCols = 'db', ...)
        }))
    }
    if(is.null(gps) ||
       nrow(gps) == 0) {
        pamWarning('No gps data found for event ', id(x))
        return(x)
    }
    diffs <- numeric(0)
    for(d in seq_along(detectors(x))) {
        withGps <- suppressWarnings(
            addGps(x@detectors[[d]], gps, thresh, keepDiff=TRUE, ...)
        )
        hasNa <- is.na(withGps$Latitude) | is.na(withGps$Longitude)
        if(any(hasNa)) {
            diffs <- c(diffs, withGps$timeDiff[hasNa])
        }
        withGps$timeDiff <- NULL
        x@detectors[[d]] <- withGps
    }
    if(length(diffs) > 0) {
        pamWarning(length(diffs), ' GPS matches in event ', id(x), ' exceeded time threshold, setting',
                   ' to NA. (Average ', round(mean(diffs/3600, na.rm=TRUE),1), ' hours apart)',
                   which=sys.nframe()-1)
    }
    x
})

#' @rdname addGps
#' @export
#'
setMethod('addGps', 'list', function(x, gps, thresh = 3600, ...) {
    lapply(x, function(y) addGps(y, gps, thresh, ...))
})

#' @rdname addGps
#' @importFrom data.table rbindlist
#' @export
#'
setMethod('addGps', 'AcousticStudy', function(x, gps=NULL, thresh = 3600, ...) {
    now <- Sys.time()
    if(is.null(gps)) {
        dbExists <- file.exists(files(x)$db)
        if(any(!dbExists)) {
            message('Could not locate database file(s) ',
                    printN(basename(files(x)$db[!dbExists]), 6),
                    ', running "updateFiles()" may help.')
        }
        gps <- rbindlist(lapply(files(x)$db, function(db) {
            if(!file.exists(db)) {
                return(NULL)
            }
            # gpsFromDb(db, extraCols=c('Speed', 'Heading', 'MagneticVariation', 'db'), ...)
            gpsFromDb(db, extraCols = 'db', ...)
        }))
        if(is.null(gps) ||
           nrow(gps) == 0) {
            pamWarning('No gps data found in any databases.')
            return(x)
        }
        # gps <- checkGpsKey(gps)
        # # setkeyv(gps, c(key(gps), 'db'))
        # events(x) <- lapply(events(x), function(y) {
        #     thisGps <- gps[gps$db %in% files(y)$db, c('UTC', 'Latitude', 'Longitude'), with=FALSE]
        #     addGps(y, thisGps, thresh, ...)
        # })
        
    } #else {
    #     gps <- checkGpsKey(gps)
    #     events(x) <- lapply(events(x), function(y) {
    #         if('db' %in% colnames(gps)) {
    #             thisGps <- gps[gps$db %in% files(y)$db, c('UTC', 'Latitude', 'Longitude'), with=FALSE]
    #         } else {
    #             thisGps <- gps
    #         }
    #         addGps(y, thisGps, thresh, ...)
    #     })
    # }
    gps <- checkGpsKey(gps)
    events(x) <- lapply(events(x), function(y) {
        if('db' %in% colnames(gps)) {
            thisGps <- gps[gps$db %in% files(y)$db, c('UTC', 'Latitude', 'Longitude'), with=FALSE]
        } else {
            thisGps <- gps
        }
        suppressWarnings(addGps(y, thisGps, thresh, ...))
    })
    gps(x) <- gps
    warns <- mget('PAMWARNING', envir = sys.frame(1), ifnotfound = NA)[['PAMWARNING']]
    if(is.data.frame(warns)) {
        warns <- warns[warns$functionName == 'addGps' &
                           warns$time > now, ]
        if(nrow(warns) > 0) {
            pamWarning(nrow(warns), ' events had GPS matches exceeding the time',
                       ' threshold, use "getWarnings(x)" for details.')
        }
    }
    x <- .addPamWarning(x)
    x
})

#' @rdname addGps
#' @export
#'
setMethod('addGps', 'ANY', function(x, gps, thresh = 3600, ...) {
    warning('No addGps method for object type', class(x))
    x
})

globalVariables(c('UTC', 'gpsTime', 'dataTime', 'db'))

#' @importFrom data.table :=
#' @importFrom PAMmisc addPgGps
#'
gpsFromDb <- function(db, extraCols=NULL, bounds=NULL) {
    con <- dbConnect(db, drv=SQLite())
    on.exit(dbDisconnect(con))
    if(!('gpsData' %in% dbListTables(con))) {
        cat('No "gpsData" table found in database', basename(db), '\n')
        addChoice <- menu(title = 'Would you like to add gps data to this database?',
                          choices = c('Yes', 'No'))
        if(addChoice %in% c(0, 2)) {
            return(NULL)
        }
        thisAddGps <- functionParser(addPgGps, skipArgs = c('db', 'gps'))
        cat(paste0('Select a GPS file for database ', basename(db)))
        chooseGps <- tk_choose.files(caption = paste0('Select a GPS file for database ', basename(db),':'))
        if(length(chooseGps) == 0) {
            return(NULL)
        }
        doAdd <- thisAddGps(db=db, gps=chooseGps)
        return(gpsFromDb(db, extraCols, bounds))
    }
    thisGps <- dbReadTable(con, 'gpsData')
    setDT(thisGps)
    # thisGps[, db := db]
    thisGps$db <- db
    if(!is.null(extraCols)) {
        extraCols <- extraCols[extraCols %in% colnames(thisGps)]
    }
    thisGps <- thisGps[, c('UTC', 'UTCMilliseconds', 'Latitude', 'Longitude', extraCols), with=FALSE]
    # thisGps[, UTC := pgDateToPosix(UTC) + UTCMilliseconds / 1e3]
    thisGps$UTC <- pgDateToPosix(thisGps$UTC) + thisGps$UTCMilliseconds / 1e3
    if(!is.null(bounds)) {
        thisGps <- thisGps[UTC <= bounds[2] & UTC >= bounds[1], ]
    }
    thisGps
}

checkGpsKey <- function(gps) {
    if(!inherits(gps, 'data.frame')) {
        stop('"gps" must be a data.frame or data.table.')
    }
    if(!inherits(gps, 'data.table')) {
        setDT(gps)
    }
    if('Channel' %in% colnames(gps)) {
        if(is.null(key(gps)) ||
           !all(key(gps) %in% c('Channel', 'UTC'))) {
            setkeyv(gps, c('Channel', 'UTC'))
        }
    } else {
        if(is.null(key(gps)) ||
           !(key(gps) %in% 'UTC')) {
            setkeyv(gps, 'UTC')
        }
    }
    gps
}
