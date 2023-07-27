#' @title Add Measures
#'
#' @description Adds "measures" to an AcousticStudy or AcousticEvent.
#'   A "measure" is an event-level variable that will be exported
#'   alongside data from that event
#'
#' @param x an \linkS4class{AcousticStudy} or
#'   \linkS4class{AcousticEvent} object
#' @param measures the measures to add. Can either be a named list,
#'   where names match event names of \code{x} or a dataframe with
#'   column \code{eventId} matching the event names of \code{x}. If
#'   a list, every item within the list must also be named by the
#'   variable name. All
#'   other data within \code{measures} will be added as new measures
#' @param replace logical flag whether or not to replace
#'
#' @author Taiki Sakai \email{taiki.sakai@@noaa.gov}
#'
#' @return object of same class as \code{x} with measures added
#'
#' @examples
#'
#' data(exStudy)
#' measList <- list('Example.OE1' = list(a=1, b=2),
#'                  'Example.OE2' = list(a=2, b=3)
#'                  )
#' exMeasure <- addMeasures(exStudy, measList)
#' print(getMeasures(exMeasure))
#' measDf <- data.frame(eventId = c('Example.OE1', 'Example.OE2'),
#'                      a=4:5,
#'                      b=6:7)
#' exMeasure <- addMeasures(exMeasure, measDf, replace=TRUE)
#' getMeasures(exMeasure)
#'
#' @export
#'
addMeasures <- function(x, measures, replace=TRUE) {
    if(!is.AcousticStudy(x) & !is.AcousticEvent(x)) {
        stop('"x" must be an AcousticStudy or AcousticEvent')
    }
    measures <- checkMeasures(x, measures)
    if(is.AcousticEvent(x)) {
        thisMeas <- measures[measures$eventId == id(x), ]
        thisMeas <- dropCols(thisMeas, 'eventId')
        oldMeas <- ancillary(x)$measures
        if(!is.null(oldMeas)) {
            thisMeas <- safeListAdd(oldMeas, thisMeas, replace=replace)
        }
        ancillary(x)$measures <- thisMeas
        return(x)
    }
    for(e in seq_along(events(x))) {
        thisMeas <- measures[measures$eventId == id(x[[e]]), ]
        if(is.null(thisMeas) || nrow(thisMeas) == 0) {
            next
        }
        thisMeas <- dropCols(thisMeas, 'eventId')
        oldMeas <- ancillary(x[[e]])$measures
        if(!is.null(oldMeas)) {
            thisMeas <- safeListAdd(oldMeas, thisMeas, replace=replace)
        }
        ancillary(x[[e]])$measures <- thisMeas
    }
    x
}

checkMeasures <- function(x, measures) {
    isEvent <- is.AcousticEvent(x)
    if(isEvent) {
        evNames <- id(x)
    } else {
        evNames <- sapply(events(x), id)
    }

    if(inherits(measures, 'list')) {
        if(!any(evNames %in% names(measures))) {
            stop('Names of "measures" did not contain any of the event names')
        }
        measures <- bind_rows(measures, .id='eventId')
    }
    if(!inherits(measures, 'data.frame')) {
        stop('"measures" must be a list or dataframe')
    }
    if(!'eventId' %in% colnames(measures) ||
       !any(evNames %in% measures$eventId)) {
        stop('"eventId" column of "measures" did not contain any of the event names')
    }
    measures
}

#' @rdname addMeasures
#' @export
#'
getMeasures <- function(x) {
    if(is.AcousticEvent(x)) {
        meas <- ancillary(x)$measures
        return(
            data.frame(eventId=id(x),
                       meas)
        )
    }
    if(!is.AcousticStudy(x)) {
        stop('"x" must be an AcousticEvent or AcousticStudy')
    }
    bind_rows(
        lapply(events(x), getMeasures),
        .id='eventId'
    )
}
