#' @title Calculate Inter-Click Interval
#'
#' @description Calculate inter-click interval for click data
#'
#' @param x a \linkS4class{AcousticStudy} object, a list of \linkS4class{AcousticEvent}
#'   objects, or a single \linkS4class{AcousticEvent} object
#' @param time the time measurement to use. \code{start} will use the \code{UTC} value,
#'   \code{peak} will use the \code{peakTime} value if present (currently present in
#'   \code{standardClickCalcs}, this is the time of the peak of the waveform)
#' @param iciRange optional range of allowed ICI (time to next detection) values.
#'   Values outside of this range will be removed before calculating the modal ICI.
#' @param callType the call type to calculate ICI for, usually this is \code{click}
#'   but also allows users to specify \code{whistle} or \code{cepstrum} to calculate this
#'   using other detector data
#' @param verbose logical flag to print messages
#' @param \dots not currently used
#'
#' @details Calculates the ICI for each individual detector and across all detectors.
#'   ICI calculation is done by ordering all individual detections by time, then taking
#'   the difference between consecutive detections and approximating the mode value.
#'
#' @return the same object as \code{x}, with ICI data added to the "ancillary" slot
#'   of each AcousticEvent. Two items will be added. $ici contains all of the
#'   individual inter-click intervals used to calculate the ICI, as well as an "All"
#'   ICI using all the combined data. $measures will also have a ICI measurement added
#'   for each detector, this will be the single modal value. Data in the $measures spot
#'   can be exported easily to modeling algorithms. \code{getICI} will just return either
#'   the values stored in $measures for \code{type = 'value'} or a dataframe of the
#'   individual ICI values used to calculate these (with columns indicating separate Channels,
#'   eventIds, and detectorNames) for \code{type = 'data'}
#'
#' @author Taiki Sakai \email{taiki.sakai@@noaa.gov}
#'
#' @examples
#'
#' # setting up example data
#' data(exStudy)
#' exStudy <- calculateICI(exStudy)
#' # each event has its ICI data stored separately, these are 0
#' # because there is only a single click in this event
#' ancillary(exStudy[[1]])$ici
#' # also saves it in measures that will get exported for modeling
#' ancillary(exStudy[[1]])$measures
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
                                    iciRange=c(0, Inf),
                                    callType = c('click', 'whistle', 'cepstrum', 'gpl'),
                                    verbose=TRUE,
                                    ...) standardGeneric('calculateICI'))

#' @rdname calculateICI
#' @export
#'
setMethod('calculateICI', 'AcousticStudy', function(x,
                                                    time=c('UTC', 'peakTime'),
                                                    iciRange=c(0, Inf),
                                                    callType = c('click', 'whistle', 'cepstrum', 'gpl'),
                                                    verbose=TRUE,
                                                    ...) {
    events(x) <- lapply(events(x), function(e) calculateICI(e, time=time, iciRange=iciRange,
                                                            callType=callType, verbose=verbose, ...))
    x
})

#' @rdname calculateICI
#' @export
#'
setMethod('calculateICI', 'AcousticEvent', function(x,
                                                    time=c('UTC', 'peakTime'),
                                                    iciRange=c(0, Inf),
                                                    callType = c('click', 'whistle', 'cepstrum', 'gpl'),
                                                    verbose=TRUE,
                                                    ...) {
    callType <- match.arg(callType)
    detData <- getDetectorData(x)[[callType]]
    time <- match.arg(time)
    detData <- detData[!is.na(detData[[time]]), ]
    if(is.null(detData)) {
        if(verbose) {
            message('No detector data found for call type "', callType, '" in event "', id(x), '"\n')
        }
        return(x)
    }
    detNames <- unique(detData$detectorName)
    iciList <- vector('list', length = length(detNames) + 1)
    names(iciList) <- c(detNames, 'All')
    for(d in detNames) {
        thisIci <- dfTimeToNext(detData[detData$detectorName == d, ], time)
        thisIci$detectorName <- d
        iciList[[d]] <- thisIci
    }
    allIci <- dfTimeToNext(detData, time)
    allIci$detectorName <- 'All'
    iciList[['All']] <- allIci

    oldIci <- ancillary(x)$ici
    if(!is.null(oldIci)) {
        iciList <- safeListAdd(oldIci, iciList)
    }
    ancillary(x)$ici <- iciList
    # browser()
    # filter outliers before mode calc
    iciMode <- lapply(iciList, function(i) {
        calcIciMode(i$ici, iciRange=iciRange)
    })
    names(iciMode) <- paste0(names(iciMode), '_ici')
    ancillary(x)$measures <- safeListAdd(ancillary(x)$measures, iciMode)
    x <- .addPamWarning(x)
    x
})

calcIciMode <- function(ici, iciRange=c(0, Inf)) {
    if(length(iciRange) != 2) {
        iciRange <- c(0, Inf)
    }
    inRange <- ici > iciRange[1] &
        ici < iciRange[2]
    ici <- ici[inRange]
    if(length(ici) == 0) {
        return(0)
    }
    if(length(ici) == 1) {
        return(ici)
    }
    if(!is.na(sd(ici)) &&
       sd(ici) != 0) {
        iciZ <- (ici - mean(ici)) / sd(ici)
        ici <- ici[abs(iciZ) < 2]
        if(any(abs(iciZ) < 2)) {
            iciZ <- (ici - mean(ici)) / sd(ici)
            ici <- ici[abs(iciZ) < 2]
        }
        if(length(ici) == 0) {
            return(0)
        }
        if(length(ici) == 1) {
            return(ici)
        }
    }
    den <- density(ici)
    mode <- den$x[which.max(den$y)]
    if(mode < 0) mode <- 0
    mode
}

dfTimeToNext <- function(x, time='UTC') {
    # check if peakTime is full time or just time within the waveform
    # 1e4 arbitrary, but UTC as numeric will be ~ 1e9
    if(time == 'peakTime' &&
       is.numeric(x$peakTime[1]) &&
       x$peakTime[1] < 1e4) {
        x$peakTime <- as.numeric(x$UTC) + x$peakTime
    }
    if('Channel' %in% colnames(x)) {
        bind_rows(lapply(unique(x$Channel), function(c) {
            calcTimeToNext(x[x$Channel == c, c('UID', 'BinaryFile', 'Channel', time)], time)
        }))
    } else {
        calcTimeToNext(x[, c('UID', 'BinaryFile',time)], time)
    }
}

calcTimeToNext <- function(x, time) {
    if(nrow(x) == 1) {
        x$ici <- 0
        return(x)
    }
    x$sort <- as.numeric(x[[time]])
    x <- arrange(x, .data$sort)
    # ici <- c(x$sort[2:nrow(x)], x$sort[nrow(x)]) - x$sort
    # x$ici <- ici
    x$ici <- c(diff(x$sort), 0)
    x$sort <- NULL
    x
}

#' @export
#' @rdname calculateICI
#' @param type the type of data to return, one of 'value' or 'data'. 'value' returns
#'  the single ICI value for each detector, 'data' returns all the individual ICI values
#'  used to calculate the number returned by 'value'
#'
getICI <- function(x, type=c('value', 'data')) {
    type <- match.arg(type)
    if(is.AcousticStudy(x)) {
        result <- suppressPamWarnings(lapply(events(x), function(e) getICI(e, type)))
        noICI <- sapply(result, function(r) is.null(r))
        if(all(noICI)) {
            pamWarning('No ICI data found, run "calculateICI" first')
            return(NULL)
        }
        if(any(noICI)) {
            pamWarning('No ICI data found in event(s) ', names(result)[noICI], n=6)
        }
        if(type == 'data') {
            result <- bind_rows(result, .id='eventId')
        }
        return(result)
    }
    if(!is.AcousticEvent(x)) {
        stop('Input must be an AcousticStudy or AcousticEvent')
    }
    switch(type,
           'value' = {
               isIci <- grep('_ici$', names(ancillary(x)$measures), value=TRUE)
               if(length(isIci) == 0) {
                   pamWarning('No ICI data found in event ', id(x), ' run "calculateICI" first.')
                   return(NULL)
               }
               ancillary(x)$measures[isIci]
           },
           'data' = {
               if(is.null(ancillary(x)$ici)) {
                   pamWarning('No ICI data found in event ', id(x), ' run "calculateICI" first.')
                   return(NULL)
               }
               bind_rows(ancillary(x)$ici)
           })
}

