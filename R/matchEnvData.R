#' @title Match Environmental Data to an AcousticStudy Object
#'
#' @description Extracts all variables from a netcdf file matching Longitude,
#'   Latitude, and UTC coordinates of the start of each AcousticEvent object.
#'   Matched values are stored in the "ancillary" slot of each event
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
#' @param progress logical flag to show progress bar
#' @param \dots other parameters to pass to \link[PAMmisc]{ncToData}
#'
#' @return original data object with environmental data added to the \code{ancillary} slot
#'   of each event. Complete data will be stored in \code{ancillary(data)$environmental},
#'   and the mean of each downloaded variable will be stored in \code{ancillary(data)$measures}
#'   so that it can be exported for modeling. For each event the coordinates associated with
#'   the earliest UTC value in that event are used to match
#'
#' @examples
#'
#' data(exStudy)
#' nc <- system.file('extdata', 'sst.nc', package='PAMmisc')
#' # suppressing warnings because nc coordinates dont align with this data,
#' # function warns of possible coordinate mismatch
#' exStudy <- suppressWarnings(matchEnvData(exStudy, nc=nc, progress=FALSE))
#' str(ancillary(exStudy[[1]])$environmental)
#' ancillary(exStudy[[1]])$measures
#'
#' @author Taiki Sakai \email{taiki.sakai@@noaa.gov}
#' @rdname matchEnvData
#' @importMethodsFrom PAMmisc matchEnvData
#' @export
#'
setMethod('matchEnvData',
          'AcousticEvent',
          function(data, nc=NULL, var=NULL, buffer=c(0,0,0), FUN = c(mean, median, sd),
                   fileName = NULL, progress=TRUE, ...) {
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
              envData <- matchEnvData(eventStart, nc=nc, var=var, buffer=buffer, FUN=FUN,
                                      fileName=fileName, progress=progress, ...)
              addEnvToEvent(data, envData)
          }
)

#' @export
#' @rdname matchEnvData
#'
setMethod('matchEnvData',
          'AcousticStudy',
          function(data, nc=NULL, var=NULL, buffer=c(0,0,0), FUN = c(mean, median, sd),
                   fileName = NULL, progress=TRUE, ...) {
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
              envData <- matchEnvData(eventStart, nc=nc, var=var, buffer=buffer, FUN=FUN,
                                      fileName=fileName, progress=progress, ...)
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
