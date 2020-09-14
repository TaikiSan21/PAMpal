#' @title Match Data From an Existing Netcdf File or Download and Match
#'
#' @description Extracts all variables from a netcdf file matching Longitude,
#'   Latitude, and UTC coordinates in given dataframe
#'
#' @param data an \linkS4class{AcousticStudy} or \linkS4class{AcousticEvent} object
#'   that must have GPS data added to it using the \link{addGps} functions
#' @param nc name of a netcdf file, ERDDAP dataset id, or an edinfo object
#' @param var (optional) vector of variable names
#' @param buffer vector of Longitude, Latitude, and Time (seconds) to buffer around
#'   each datapoint. All values within the buffer will be used to report the mean,
#'   median, and standard deviation
#' @param FUN a vector or list of functions to apply to the data. Default is to apply
#'   mean, median, and standard deviation calculations
#' @param fileName (optional) file name to save downloaded nc file to. If not provided,
#'   then no nc files will be stored, instead small temporary files will be downloaded
#'   and then deleted. This can be much faster, but means that the data will need to be
#'   downloaded again in the future. If \code{fileName} is provided, then the function
#'   will attempt to download a single nc file covering the entire range of your data.
#'   If your data spans a large amount of time and space this can be problematic.
#' @param \dots other parameters to pass to \link[PAMmisc]{ncToData}
#'
#' @return original data object with environmental data added to the \code{ancillary} slot
#'   of each event. Complete data will be stored in \code{ancillary(data)$environmental},
#'   and the mean of each downloaded variable will be stored in \code{ancillary(data)$measures}
#'   so that it can be exported for modeling
#'
#' @author Taiki Sakai \email{taiki.sakai@@noaa.gov}
#' @rdname matchEnvData
#' @importMethodsFrom PAMmisc matchEnvData
#' @export
#'
setMethod('matchEnvData',
          'AcousticEvent',
          function(data, nc=NULL, var=NULL, buffer=c(0,0,0), FUN = c(mean, median, sd), fileName = NULL, ...) {
              if(is.list(FUN) &&
                 is.null(names(FUN))) {
                  names(FUN) <- as.character(substitute(FUN))[-1]
              } else if(is.function(FUN)) {
                  tmpName <- as.character(substitute(FUN))
                  FUN <- list(FUN)
                  names(FUN) <- tmpName
              }
              eventStart <- getEventStart(data)
              if(is.null(eventStart)) {
                  return(data)
              }
              envData <- matchEnvData(eventStart, nc=nc, var=var, buffer=buffer, FUN=FUN, fileName=fileName, ...)
              addEnvToEvent(data, envData)
          }
)

#' @export
#' @rdname matchEnvData
#'
setMethod('matchEnvData',
          'AcousticStudy',
          function(data, nc=NULL, var=NULL, buffer=c(0,0,0), FUN = c(mean, median, sd), fileName = NULL, ...) {
              eventStart <- do.call(rbind, lapply(events(data), function(x) {
                  getEventStart(x)
              }))
              if(is.null(eventStart)) {
                  return(data)
              }
              if(is.list(FUN) &&
                 is.null(names(FUN))) {
                  names(FUN) <- as.character(substitute(FUN))[-1]
              } else if(is.function(FUN)) {
                  tmpName <- as.character(substitute(FUN))
                  FUN <- list(FUN)
                  names(FUN) <- tmpName
              }
              envData <- matchEnvData(eventStart, nc=nc, var=var, buffer=buffer, FUN=FUN, fileName=fileName, ...)
              for(e in seq_along(events(data))) {
                  events(data)[[e]] <- addEnvToEvent(events(data)[[e]], envData)
              }
              data
          }
)

addEnvToEvent <- function(event, env) {
    env <- env[env$event == id(event), ]
    envList <- as.list(env)
    oldEnv <- ancillary(event)$environmental
    if(!is.null(oldEnv)) {
        envList <- safeListAdd(oldEnv, envList)
    }
    ancillary(event)$environmental <- envList
    measureOnly <- !grepl('_median$|_sd$|Longitude|match|Latitude|UTC|event', names(env))
    measList <- env[measureOnly]
    oldMeas <- ancillary(event)$measures
    if(!is.null(oldMeas)) {
        measList <- safeListAdd(oldMeas, measList)
    }
    ancillary(event)$measures <- measList
    event
}

getEventStart <- function(data) {
    dets <- detectors(data)
    if(length(dets) == 0) {
        return(NULL)
    }
    hasDets <- sapply(dets, function(x) nrow(x) > 0)
    if(sum(hasDets) == 0) {
        return(NULL)
    }
    dets <- dets[hasDets]
    if(!all(c('UTC', 'Longitude', 'Latitude') %in% colnames(dets[[1]]))) {
        warning(paste0('Event ', id(data), ' does not have GPS data added.'))
        return(NULL)
    }
    coordsOnly <- do.call(rbind, lapply(dets, function(x) {
        x[, c('UTC', 'Longitude', 'Latitude')]
    }))
    coordsOnly$event <- id(data)
    # for now use earliest time in event, first row only just in case a tie
    coordsOnly[which.min(coordsOnly$UTC), ][1, ]
}
