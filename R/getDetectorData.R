#' @title Extract and Combine Detector Data
#'
#' @description Extracts just the detector data from all of \code{x}, and will
#'   combine all detections from each call type (currently whistle, click,
#'   and cepstrum) into a single data frame.
#'
#' @param x data to extract detector data from, either an \code{AcousticStudy},
#'   \code{AcousticEvent} or list of \code{AcousticEvent} object
#'
#' @details The purpose of this function is to extract your data out of
#'   \code{PAMpal}'s S4 classes and put them into an easier format to work with.
#'   The output will be a list of up to three data frames, one for each call type
#'   found in your data. Each different call type will have had different processing
#'   applied to it by \code{processPgDetections}. Additionally, each detector will
#'   have its associated event id, the name of the detector, and the species id
#'   attached to it (species will be \code{NA} if not set). All detections from each
#'   call type will be combined into a single large data frame
#'
#' @return A list of data frames containing all detection data from \code{x},
#'   named by call type ('click', 'whistle', or 'cepstrum').
#'
#' @author Taiki Sakai \email{taiki.sakai@@noaa.gov}
#'
#' @examples
#'
#' data(exStudy)
#' dets <- getDetectorData(exStudy)
#' names(dets)
#' str(dets$click)
#' # works on single events as well
#' oneDets <- getDetectorData(exStudy[[1]])
#' str(oneDets$click)
#'
#' @importFrom PAMmisc squishList
#' @export
#'
getDetectorData <- function(x) {
    if(is.data.frame(x)) {
        return(x)
    }
    if(is.AcousticStudy(x)) {
        return(getDetectorData(events(x)))
    }
    if(is.list(x)) {
        if(length(x) == 0) {
            return(NULL)
        }
        result <- lapply(x[sapply(x, is.AcousticEvent)], function(e) {
            getDetectorData(e)
        })
        names(result) <- NULL
        result <- unlist(result, recursive=FALSE)
        return(squishList(result))
    }
    # base case one acev
    dets <- detectors(x)
    callTypes <- getCallType(dets)
    for(d in seq_along(dets)) {
        dets[[d]]$eventId <- id(x)
        dets[[d]]$detectorName <- names(dets)[d]

        if(is.null(species(x)$id)) {
            dets[[d]]$species <- NA_character_
        } else {
            dets[[d]]$species <- species(x)$id[1]
        }
    }
    names(dets) <- callTypes
    squishList(dets)
}

# possible call types are whistles clicks cepstrum
getCallType <- function(x) {
    if(is.data.frame(x)) {
        return(attr(x, 'calltype'))
    }
    if(is.AcousticEvent(x)) {
        x <- detectors(x)
    }
    if(is.list(x) &&
       all(sapply(x, is.data.frame))) {
        return(
            sapply(x, function(d) {
                attr(d, 'calltype')
            })
        )
    }
    NULL
}

#' @export
#' @rdname getDetectorData
#'
getClickData <- function(x) {
    getDetectorData(x)$click
}

#' @export
#' @rdname getDetectorData
#'
getWhistleData <- function(x) {
    getDetectorData(x)$whistle
}

#' @export
#' @rdname getDetectorData
#'
getCepstrumData <- function(x) {
    getDetectorData(x)$cepstrum
}

