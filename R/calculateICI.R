#' @title Calculate Inter-Click Interval
#'
#' @description Calculate inter-click interval for click data
#'
#' @param x a \linkS4class{AcousticStudy} object, a list of \linkS4class{AcousticEvent}
#'   objects, or a single \linkS4class{AcousticEvent} object
#' @param time the time measurement to use. \code{start} will use the \code{UTC} value,
#'   \code{peak} will use the \code{peakTime} value if present (currently present in
#'   \code{standardClickCalcs}, this is the time of the peak of the waveform)
#' @param callType the call type to calculate ICI for, usually this is \code{click}
#'   but also allows users to specify \code{whistle} or \code{cepstrum} to calculate this
#'   using other detector data
#' @param \dots not used
#'
#' @details Calculates the ICI for each individual detector and across all detectors.
#'   ICI calculation is done by ordering all individual detections by time, then taking
#'   the difference between consecutive detections and taking the mode value.
#'
#' @return ICI data
#'
#' @author Taiki Sakai \email{taiki.sakai@@noaa.gov}
#'
#' @name calculateICI
#'
#' @importFrom dplyr bind_rows
#' @importFrom stats density
#'
#' @export
#'
setGeneric('calculateICI', function(x,
                                    time=c('UTC', 'peakTime'),
                                    callType = c('click', 'whistle', 'cepstrum'),
                                    ...) standardGeneric('calculateICI'))


#' @rdname calculateICI
#' @export
#'
setMethod('calculateICI', 'AcousticStudy', function(x,
                                                    time=c('UTC', 'peakTime'),
                                                    callType = c('click', 'whistle', 'cepstrum'),
                                                    ...) {
  events(x) <- lapply(events(x), function(e) calculateICI(e, time, callType, ...))
  x
})

#' @rdname calculateICI
#' @export
#'
setMethod('calculateICI', 'AcousticEvent', function(x,
                                                    time=c('UTC', 'peakTime'),
                                                    callType = c('click', 'whistle', 'cepstrum'),
                                                    ...) {
  callType <- match.arg(callType)
  detData <- getDetectorData(x)[[callType]]
  if(is.null(detData)) {
    cat('No detector data found for call type', callType, 'in event', id(x), '\n')
    return(x)
  }
  detNames <- unique(detData$detectorName)
  iciList <- vector('list', length = length(detNames) + 1)
  names(iciList) <- c(detNames, 'All')
  time <- match.arg(time)
  for(d in detNames) {
      iciList[[d]] <- data.frame(name=d,
                                 ici=dfICI(detData[detData$detectorName == d, ], time),
                                 stringsAsFactors = FALSE)
  }
  iciList[['All']] <- data.frame(name='All',
                                 ici=dfICI(detData, time),
                                 stringsAsFactors = FALSE)
  oldIci <- ancillary(x)$ici
  if(!is.null(oldIci)) {
      iciList <- safeListAdd(oldIci, iciList)
  }
  ancillary(x)$ici <- iciList
  # browser()
  iciMode <- lapply(iciList, function(i) {
      den <- density(i$ici)
      mode <- den$x[which.max(den$y)]
      if(mode < 0) mode <- 0
      mode
  })
  names(iciMode) <- paste0(names(iciMode), '_ici')
  ancillary(x)$measures <- safeListAdd(ancillary(x)$measures, iciMode)
  x
})

dfICI <- function(x, time='UTC', plot=FALSE) {
    # check if peakTime is full time or just time within the waveform
    # 1e4 arbitrary, but UTC as numeric will be ~ 1e9
    if(time == 'peakTime' &&
       is.numeric(x$peakTime[1]) &&
       x$peakTime[1] < 1e4) {
        x$peakTime <- as.numeric(x$UTC) + x$peakTime
    }
    if('Channel' %in% colnames(x)) {
        unlist(lapply(unique(x$Channel), function(c) {
            calcICI(x[x$Channel == c, time])
        }))
    } else {
        calcICI(x[, time])
    }
}


calcICI <- function(x) {
    if(length(x) == 1) return(0)
    time <- sort(as.numeric(x))
    ici <- time - c(time[1], time[1:(length(time)-1)])
    ici
}

