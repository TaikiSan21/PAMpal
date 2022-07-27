#' @title Add Hydrophone Depth Data to an AcousticStudy
#'
#' @description Add hydrophone depth to an AcousticStudy or AcousticEvent
#'
#' @param x an \linkS4class{AcousticStudy} to add depth data to
#' @param depth a data frame of depth values to match to data from \code{x}.
#'   Must have column \code{UTC}, and a column containing depth data to be
#'   specified by \code{depthCol}. If not provided and \code{x} is an
#'   \linkS4class{AcousticEvent} or \linkS4class{AcousticStudy} object, then
#'   the depth data will be read from the databases contained in the \code{files}
#'   slot of \code{x}
#' @param depthCol the name of the column containing depth in the dataframe or
#'   database. If left as \code{NULL}, will search for a single column containing
#'   the word "depth" or "Depth"
#' @param thresh maximum time apart in seconds for matching depth to
#'   data, if the closest value is more than \code{thresh} apart then the
#'   depth value will be set to \code{NA}
#' @param \dots additional arguments for other methods
#'
#' @details Depth values will be matched to the data
#'   by using data.table's rolling join with \code{roll='nearest'}. After the
#'   join is done, the time difference between the matched rows is checked
#'   and any that are greater than the set threshold are set to NA. This is
#'   done to prevent accidentally matching weird things if an incomplete set
#'   of depth data is provided.
#'
#'   If \code{x} is an \linkS4class{AcousticEvent} or \linkS4class{AcousticStudy},
#'   then \code{depth} can be omitted and will be read from the databases contained
#'   in the \code{files} slot of \code{x}.
#'
#' @return the same data as \code{x}, with depth data added. All AcousticEvents will
#'   have depth data added to all detector dataframes as column \code{hpDepth}
#'
#' @author Taiki Sakai \email{taiki.sakai@@noaa.gov}
#'
#' @examples
#'
#' data(exStudy)
#' # need to update database file to local directory
#' db <- system.file('extdata', 'Example.sqlite3', package='PAMpal')
#' exStudy <- updateFiles(exStudy, db=db, bin=NA, verbose=FALSE)
#' exStudy <- addHydrophoneDepth(exStudy)
#' getClickData(exStudy[1])
#'
#' @importFrom data.table data.table setkeyv key setDT setDF as.data.table
#' @export
#'
addHydrophoneDepth <- function(x, depth=NULL, depthCol=NULL, thresh=60, ...) {
  if(!is.null(depth) &&
     !is.data.frame(depth)) {
    stop('"depth" must be a dataframe or NULL to read from database.')
  }
  if(is.data.frame(depth)) {
    needCols <- c('UTC', depthCol)
    missingCols <- needCols[!(needCols %in% colnames(depth))]
    if(length(missingCols) > 0) {
      stop('Depth data needs column(s) named ', paste0(missingCols, collapse=', '))
    }
    if(!inherits(depth$UTC, 'POSIXct')) {
      stop('UTC must be converted to POSIXct in "depth" data')
    }
    depthCol <- checkDepthCol(depth, depthCol)
    depth <- select(depth, all_of(c('UTC', depthCol)))
    for(e in seq_along(events(x))) {
      x[[e]] <- depthToEvent(x[[e]], depth, thresh)
    }
    x <- .addPamWarning(x)
    return(x)
  }
  dbMap <- bind_rows(lapply(events(x), function(e) {
    list(db=files(e)$db,
         id=id(e))
  }))
  for(d in unique(dbMap$db)) {
    thisIds <- dbMap$id[dbMap$db == d]
    thisDepth <- depthFromDb(d, depthCol)
    for(i in thisIds) {
      x[[i]] <- depthToEvent(x[[i]], thisDepth, thresh)
    }
  }
  x <- .addPamWarning(x)
  x
}

depthToEvent <- function(event, depth, thresh) {
  for(d in seq_along(detectors(event))) {
    detectors(event)[[d]] <- depthToDf(detectors(event)[[d]], depth, thresh)
  }
  event
}

depthToDf <- function(x, depth, thresh) {
  # Check for right columns and proper types
  needCols <- c('UTC', 'depth')
  missingCols <- needCols[!(needCols %in% colnames(depth))]
  thisType <- attr(x, 'calltype')
  if(length(missingCols) > 0) {
    pamWarning('Depth data needs column(s) named ', paste(missingCols, collapse = ', '))
    return(x)
  }
  if(!('UTC' %in% colnames(x))) {
    pamWarning('Data needs column UTC.')
    return(x)
  }
  depth <- rename(depth, 'hpDepth' = depth)
  if(!('POSIXct' %in% class(x$UTC))) x$UTC <- pgDateToPosix(x$UTC)
  # dummies for calculating time difference for threshold check later

  # setDT(x)
  x <- dropCols(x, 'hpDepth')
  x <- as.data.table(x)
  x$dataTime <- x$UTC
  # depth <- checkGpsKey(depth)
  if(!inherits(depth, 'data.table')) {
    setDT(depth)
  }
  depth$depthTime <- depth$UTC

  setkeyv(x, 'UTC')
  setkeyv(depth, 'UTC') # removing channel key from gps if its there i guess

  result <- depth[x, roll='nearest']
  result[abs(dataTime - depthTime) > thresh, c('hpDepth') := NA]
  result$UTC <- result$dataTime
  result[, c('depthTime', 'dataTime') := NULL]
  if(any(is.na(result$hpDepth))) {
    pamWarning('Some depth matches exceeded time threshold, setting',
               'value to NA.')
  }
  attr(result, 'calltype') <- thisType
  # setDF(x)
  setDF(result)
  result
}

depthFromDb <- function(db, depthCol=c('Sensor_0_Depth'), extraCols=NULL) {
  con <- dbConnect(db, drv=SQLite())
  on.exit(dbDisconnect(con))
  if(!('Hydrophone_Depth_Data' %in% dbListTables(con))) {
    pamWarning('No "Hyrophone_Depth_Data" table found indatabase', basename(db))
    return(NULL)
  }
  thisDepth <- dbReadTable(con, 'Hydrophone_Depth_Data')
  depthCol <- checkDepthCol(thisDepth, depthCol)
  thisDepth <- select(thisDepth, any_of(c('UTC', depthCol, extraCols)))
  # setDT(thisDepth)
  # thisDepth$db <- db
  thisDepth$UTC <- pgDateToPosix(thisDepth$UTC)
  thisDepth <- rename(thisDepth, 'depth' = depthCol)
  thisDepth
}

globalVariables(c('depthTime'))

checkDepthCol <- function(df, depthCol=NULL) {
  if(is.null(depthCol) ||
     !(depthCol %in% colnames(df))) {
    maybeDepth <- grep('[Dd]epth', colnames(df), value=TRUE)
    if(length(maybeDepth) == 0) {
      pamWarning('Could not find any depth column in table')
      return(NULL)
    }
    if(length(maybeDepth) > 1) {
      naCol <- apply(df[maybeDepth], 2, function(x) all(is.na(x)))
      if(sum(!naCol) == 1) {
        depthCol <- maybeDepth[!naCol]
      } else {
        pamWarning('Could not determine single depth column, please specify with "depthCol"')
        return(NULL)
      }
    }
    if(length(maybeDepth) == 1) {
      depthCol <- maybeDepth
    }
  }
  depthCol
}