#' @title Add FPOD Detector to an AcousticStudy
#'
#' @description Adds data from FPOD detector CSV files to an
#'   \linkS4class{AcousticStudy} object as new detectors of type
#'   \code{"fpod"}
#'
#' @param x an \linkS4class{AcousticStudy} object
#' @param fpod path(s) to CSV files containing FPOD detector output
#' @param detectorName name for the detector, the default \code{'FPOD'} should be
#'   fine unless you want to differentiate between multiple FPOD detectors
#'
#' @details FPOD detections are added to events based on their times. All detections
#'   between the start and end times events. NOTE: most \code{PAMpal} functions were designed
#'   with only PAMGuard data in mind, there is a chance that adding FPOD detections will
#'   cause other advanced functionality to not work.
#'
#'   Behavior is slightly different depending
#'   on how the original AcousticStudy was created. For those processed with
#'   \code{mode='db'}, the start and end times for each event are just determined
#'   by the times of detections within the event.
#'
#'   For those processed with
#'   \code{mode='recording'} or \code{mode='time'}, the start and end times for
#'   each event are determined by the start/end times of the recording files or
#'   the grouping file provided initially. This means that it is possible that
#'   there are events which initially had zero PAMGuard detections that now have
#'   FPOD detections. In these cases a new AcousticEvent will be created that only
#'   has FPOD detections, these events may not work with a variety of other \code{PAMpal}
#'   functions.
#'
#' @return the same object as \code{x} with FPOD detector data added
#'
#' @author Taiki Sakai \email{taiki.sakai@@noaa.gov}
#'
#' @export
#'
addFPOD <- function(x, fpod, detectorName='FPOD') {
    if(is.character(fpod)) {
        fpod <- loadFPOD(fpod)
    }
    grouping <- getEventGrouping(x)
    newEvent <- character(0)
    for(i in seq_along(grouping$id)) {
        # thisFpod <- filter(fpod, UTC >= grouping$start[i], UTC < grouping$end[i])
        thisFpod <- fpod[withinLHS(fpod$UTC, grouping$interval[i]), ]
        if(is.null(thisFpod) ||
           nrow(thisFpod) == 0) {
            next
        }
        if(!grouping$id[i] %in% names(events(x))) {
            newEvent <- c(newEvent, grouping$id[i])
            x[[grouping$id[i]]] <- AcousticEvent(id=grouping$id[i])
        }
        x[[grouping$id[i]]] <- addDetector(x[[grouping$id[i]]], thisFpod, name=detectorName, calltype='fpod', replace=FALSE)
    }
    if(length(newEvent) > 0) {
        pamWarning(length(newEvent), ' new events were created containing only FPOD data (',
                            'event IDs ', printN(newEvent, 6), ')')
    }
    x <- .addPamWarning(x)
    x
}

loadFPOD <- function(x) {
    fileDNE <- !file.exists(x)
    if(all(fileDNE)) {
        stop('All FPOD files provided do not exist')
    }
    if(any(fileDNE)) {
        warning(sum(fileDNE), ' FPOD files do not exist')
    }
    all <- bind_rows(lapply(x[!fileDNE], function(file) {
        fpod <- read.csv(file, stringsAsFactors = FALSE)
        if(is.null(fpod) || nrow(fpod) == 0) {
            return(NULL)
        }
        fpod$UTC <- as.POSIXct(fpod[['Time.Mn']], format='%d/%m/%Y %H:%M', tz='UTC')
        fpod$UTC <- fpod$UTC + fpod[['MuSec']] / 1e6
        colnames(fpod)[1] <- 'BinaryFile'
        toDrop <- c('Time.Mn', 'Mn')
        fpod <- dropCols(fpod, toDrop)
        # fpod$BinaryFile <- 'NA'
        fpod$UID <- paste0(format(fpod$UTC, format='%m%d%H%M'), fpod[['MuSec']])
        fpod
    }))
    all
}

getEventGrouping <- function(x) {
    if(!is.null(ancillary(x)$grouping)) {
        return(ancillary(x)$grouping)
    }
    grouping <- getTimeRange(x, mode='event', sample=FALSE)
    grouping <- bind_rows(grouping, .id='id')
    grouping$interval <- interval(grouping$start, grouping$end)
    grouping$id <- as.character(grouping$id)
    grouping
}

addDetector <- function(event, det, name=NULL, calltype=NULL, replace=FALSE) {
    if(is.null(det) ||
       nrow(det) == 0) {
        return(event)
    }
    if(is.null(name) ||
       is.null(calltype)) {
        warning('Cannot add detector without name and calltype')
        return(event)
    }
    if(isTRUE(replace)) {
        newDet <- det
    } else {
        newDet <- bind_rows(
            event@detectors[[name]],
            det
        )
    }
    attr(newDet, 'calltype') <- calltype
    event@detectors[[name]] <- newDet
    event
}
