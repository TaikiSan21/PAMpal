#' @title Add GPS Locations to Data
#'
#' @description Add GPS Lat / Long to a variety of types of data.
#'
#' @param x data to add GPS coordinates to. Must have a column \code{UTC}, and
#'   can also have an optional column \code{Channel}
#' @param gps a data frame of GPS coordinates to match to data from \code{x}.
#'   Must have columns \code{UTC}, \code{Latitude}, \code{Longitude}, and
#'   optionally \code{Channel}. If not provided and \code{x} is an
#'   \linkS4class{AcousticEvent} or \linkS4class{AcousticStudy} object, then
#'   the gps data will be read from the databases contained in the \code{files}
#'   slot of \code{x}rm(may
#' @param thresh maximum time in seconds for matching GPS coordinates to data.
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
#' @return the same data as \code{x}, with Lat/Long data added
#'
#' @author Taiki Sakai \email{taiki.sakai@@noaa.gov}
#'
#' @name addGps
#' @export
#'
setGeneric('addGps', function(x, gps=NULL, thresh = 3600, ...) standardGeneric('addGps'))

#' @rdname addGps
#' @importFrom data.table data.table setkeyv key setDT setDF
#' @importFrom dplyr mutate select
#' @importFrom magrittr %>%
#' @importFrom utils globalVariables
#' @export
#'
setMethod('addGps', 'data.frame', function(x, gps, thresh = 3600, ...) {
    # Check for right columns and proper types
    needCols <- c('UTC', 'Latitude', 'Longitude')
    missingCols <- needCols[!(needCols %in% colnames(gps))]
    if(length(missingCols) > 0) {
        warning('Gps data needs column(s) named ', paste(missingCols, collapse = ', '))
        return(x)
    }
    if(!('UTC' %in% colnames(x))) {
        warning('Data needs column UTC.')
        return(x)
    }
    if(!('POSIXct' %in% class(x$UTC))) x$UTC <- pgDateToPosix(x$UTC)
    if(!('POSIXct' %in% class(gps$UTC))) gps$UTC <- pgDateToPosix(gps$UTC)
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
    # result <- gps[x, roll='nearest'] %>%
    #     mutate(diff = abs(dataTime - gpsTime),
    #            Latitude = ifelse(diff > thresh, NA, Latitude),
    #            Longitude = ifelse(diff > thresh, NA, Longitude),
    #            UTC = dataTime) %>%
    #     select(-diff, -gpsTime, -dataTime)
    # browser()
    result <- gps[x, roll='nearest']
    result[abs(dataTime - gpsTime) > thresh, c('Latitude', 'Longitude') := NA]
    # result[, UTC := dataTime]
    result$UTC <- result$dataTime
    result[, c('gpsTime', 'dataTime') := NULL]
    if(any(is.na(result$Longitude))) {
        warning('Some GPS coordinate matches exceeded time threshold, setting',
                'value to NA.')
    }
    attr(result, 'calltype') <- thisType
    setDF(result)
})

#' @rdname addGps
#' @export
#'
setMethod('addGps', signature(x='AcousticEvent'), function(x, gps=NULL, thresh = 3600, ...) {
    if(is.null(gps)) {
        gps <- rbindlist(lapply(files(x)$db, function(db) {
            gpsFromDb(db, extraCols=c('Speed', 'Heading', 'MagneticVariation', 'db'), ...)
        }))
    }
    if(is.null(gps) ||
       nrow(gps) == 0) {
        warning('No gps data found for event ', id(x))
        return(x)
    }
    x@detectors <- addGps(x@detectors, gps, thresh, ...)
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
    if(is.null(gps)) {
        gps <- rbindlist(lapply(files(x)$db, function(db) {
            gpsFromDb(db, extraCols=c('Speed', 'Heading', 'MagneticVariation', 'db'), ...)
        }))
        if(is.null(gps) ||
           nrow(gps) == 0) {
            warning('No gps data found in any databases.')
            return(x)
        }
        gps <- checkGpsKey(gps)
        # setkeyv(gps, c(key(gps), 'db'))
        events(x) <- lapply(events(x), function(y) {
            thisGps <- gps[db %in% files(y)$db, c('UTC', 'Latitude', 'Longitude'), with=FALSE]
            addGps(y, thisGps, thresh, ...)
        })

    } else {
        gps <- checkGpsKey(gps)
        events(x) <- lapply(events(x), function(y) addGps(y, gps, thresh, ...))
    }
    gps(x) <- gps
    x
})

#' @rdname addGps
#' @export
#'
setMethod('addGps', 'ANY', function(x, gps, thresh = 3600, ...) {
    cat('No addGps method for object type', class(x))
    x
})

globalVariables(c('UTC', 'gpsTime', 'dataTime', 'db'))

#' @importFrom data.table :=
#'
gpsFromDb <- function(db, extraCols=NULL, bounds=NULL) {
    con <- dbConnect(db, drv=SQLite())
    on.exit(dbDisconnect(con))
    if(!('gpsData' %in% dbListTables(con))) {
        cat('No "gpsData" table found in database', basename(db), '\n')
        return(NULL)
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
