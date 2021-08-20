#' @title Check an AcousticStudy Object for Issues
#'
#' @description Checks for any possible issues in an \linkS4class{AcousticStudy}
#'   object, issuing warnings and saving the messages
#'
#' @details This function is called at the end of \link{processPgDetections} with
#'   default parameters, but can also be called later to investigate issues
#'   specific to each user's data. For example, if you are expecting to process data
#'   where all recordings were duty cycled to record 2 out of every 10 minutes, then
#'   setting \code{maxLength = 60*2} will alert you to any events that are longer than
#'   the 2 minute duty cycle.
#'   For continuously recorded data, the \code{maxSep} argument can be used to
#'   identify situations where there are large gaps between detections in a single
#'   event, since this could mean that detections were accidentally added to the
#'   incorrect event number during processing.
#'
#' @param x an \linkS4class{AcousticStudy} object
#' @param maxLength events with length greater than this value in seconds
#'   will trigger a warning
#' @param maxSep events containing consecutive detections greater than
#'   \code{maxSep} seconds apart will trigger a warning. This is used to
#'   check for situations where detections were possibly added to the
#'   incorrect event.
#'
#' @return returns a list of warning messages
#'
#' @examples
#'
#' data(exStudy)
#'
#' # checks if any peak frequencies are 0, so we'll force this
#' exStudy[[1]][[1]]$peak <- 0
#' checkStudy(exStudy)
#' checkStudy(exStudy, maxLength = 1, maxSep = 1)
#'
#' @export
#'
#' @author Taiki Sakai \email{taiki.sakai@@noaa.gov}
#'
checkStudy <- function(x, maxLength=Inf, maxSep=60*60*2) {
    if(length(events(x)) == 0) {
        pamWarning('No events in AcousticStudy')
        return(NULL)
    }
    dets <- getDetectorData(x)
    peak0Msg <- doCheck(checkPeakZero, dets)
    timeMsg <- doCheck(checkTime, x, length=maxLength, between=maxSep)
    NAmsg <- doCheck(checkNAVals, dets)
    list(peak0Check = peak0Msg,
         timeCheck = timeMsg,
         NACheck = NAmsg)
}

doCheck <- function(fun, x, ...) {
    msg <- 'No problems'
    tryCatch({
        msg <- fun(x=x, ...)
    },
    error = function(e) {
        return(paste0('Error: ', e$message))
    })
    msg
}
# Check for NA vals
# threshold max distance between consec detections 2 hours
#
checkNAVals <- function(x) {
    nonParamCols <- c('UID', 'UTC', 'Latitude', 'Longitude', 'BinaryFile', 'eventId',
                      'detectorName', 'species', 'Channel', 'angle', 'angleError')
    msg <- character(0)
    for(t in names(x)) {
        paramCols <- colnames(x[[t]])[!(colnames(x[[t]]) %in% nonParamCols)]
        if(is.null(x[[t]]) ||
           length(paramCols) == 0) {
            naUID <- character(0)
        } else {
            naUID <- x[[t]]$UID[apply(x[[t]][, paramCols], 1, function(i) any(is.na(i)))]
        }
        if(length(naUID) == 0) {
            next
        }
        msg <- paste0(msg, 'Found ', length(naUID), ' ', t, ' detections with NA values.\n')
    }
    if(nchar(msg) > 0) {
        pamWarning(msg)
    } else {
        msg <- 'No detections with NA values.'
    }
    msg
}

# if peak freq is 0 probably something went wrong
checkPeakZero <- function(x) {
    peak0Msg <- 'All peaks greater than 0.'
    # det <- getDetectorData(x)
    if(nrow(x$click) > 0 &&
       'peak' %in% colnames(x$click)) {
        if(any(x$click$peak == 0, na.rm=TRUE)) {
            peak0Msg <- paste0('Some clicks had a peak frequency of 0 Hz,',
                               ' consider adjusting the filter parameter',
                               ' or adding a calibration function.')
            pamWarning(peak0Msg)
        }
    }
    peak0Msg
}
# check event times are too long or have strangely long distance between
# consecutive detctions - prob means an event was mismarked
checkTime <- function(x, length=Inf, between=60*60*2) {
    evTimes <- bind_rows(lapply(events(x), function(e) {
        allTimes <- unlist(lapply(detectors(e), function(d) {
            as.numeric(d[['UTC']])
        }))
        allTimes <- sort(allTimes)
        if(length(allTimes) > 1) {
            list(evLen = diff(range(allTimes)),
                 evBtwn = max(diff(allTimes)),
                 id = id(e))
        } else {
            list(evLen = 0,
                 evBtwn = 0,
                 id = id(e)
            )
        }
    }))
    # browser()
    isLong <- evTimes$evLen > length
    isBad <- evTimes$evBtwn > between
    if(!any(isLong) &&
       !any(isBad)) {
        return('No event length issues.')
    }
    longMsg <- ''
    if(any(isLong)) {
        longIds <- paste0(evTimes$id[isLong], collapse=', ')
        longMsg <- paste0('Found ', sum(isLong), ' events longer than ',
                          length, ' seconds: ', longIds)
        pamWarning(longMsg)
    }
    badMsg <- ''
    if(any(isBad)) {
        badIds <- paste0(evTimes$id[isBad], collapse=', ')
        badMsg <- paste0('Found ', sum(isBad), ' events with detections',
                         ' more than ', between, ' seconds apart: ', badIds)
        pamWarning(badMsg)
    }
    paste0(longMsg, '\n', badMsg)
}
