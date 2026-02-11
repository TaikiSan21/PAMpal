#' @title Filter Candidate Echo Depths
#'
#' @description Filter out possible echo depths from \link{calculateEchoDepth}
#'   based on maximum depth, autocorrelation magnitude, and maximum swim
#'   speed criteria. Requires that \link{calculateEchoDepth} has been run
#'   first. This function adds a \code{keepClick} column to the data to
#'   track which detections should be used for further depth analysis by
#'   marking them as \code{FALSE} to be excluded or \code{TRUE} to be used
#'
#' @param x an \linkS4class{AcousticStudy} object that has been processed with
#'   \link{calculateEchoDepth}
#' @param time maximum time apart (seconds) for detections. Detections with no
#'   no other detection within \code{time} seconds will be marked as \code{FALSE}
#' @param depth maximum depth difference (meters) between consecutive clicks,
#'   this value should be determined by maximum swim speed
#' @param speed as an alternative to providing \code{depth}, the swim speed
#'   (meters / second) can be provided and then \code{depth} will be calculated
#'   as \code{time} * \code{speed}
#' @param maxDepth calculated depth values greater than this will be marked
#'   \code{FALSE}
#' @param minCorr detections with autocorrelation magnitude less than this
#'   will be marked as \code{FALSE}
#'
#' @author Taiki Sakai \email{taiki.sakai@@noaa.gov}
#'
#' @return the AcousticStudy \code{x} with detections marked with column
#'   \code{keepClick} as \code{TRUE} or \code{FALSE} depending if they
#'   pass the filter parameters
#'
#' @examples
#' # example not run because \link{calculateEchoDepth} must be run first,
#' # and it requires a large amount of data not stored in the package
#' \dontrun{
#' study <- calculateEchoDepth(study, wav='path/to/wavFiles')
#' study <- filterEchoDepths(study, time=30, speed=50/30, maxDepth=4000)
#' }
#'
#' @export
#'
filterEchoDepths <- function(x, time=30, depth=NULL, speed=NULL,
                         maxDepth=4000, minCorr=.01) {
    if(is.null(depth) && is.null(speed)) {
        stop('Either depth or swim speed (m/s) must be provided')
    }
    if(is.null(depth)) {
        depth <- time * speed
    }
    if(is.AcousticStudy(x) || is.AcousticEvent(x)) {
        clicks <- getClickData(x, measures=FALSE)
        detNameHolder <- select(clicks, any_of(c('UID', 'eventId', 'Channel', 'detectorName')))
        clicks <- distinct(dropCols(clicks, c('detectorName', 'ici')))
    }
    if(is.character(x)) {
        if(!file.exists(x)) {
            stop('File ', x, ' does not exist')
        }
        clicks <- read.csv(x, stringsAsFactors = FALSE)
    }
    if(is.data.frame(x)) {
        clicks <- x
    }
    clicks$UTC <- parseUTC(clicks$UTC)
    # intialize tracking column, dont keep NA depths
    if(!'keepClick' %in% colnames(clicks)) {
        clicks$keepClick <- !is.na(clicks$maxDepth)
    }
    clicks <- split(clicks, clicks$eventId)
    clicks <- lapply(clicks, function(x) {
        filterBySwimSpeed(x, time=time, depth=depth)
    })
    clicks <- bind_rows(clicks)
    clicks$keepClick[clicks[['maxDepth']] > maxDepth] <- FALSE
    clicks$keepClick[clicks[['maxMag']] < minCorr] <- FALSE
    if(is.AcousticStudy(x)) {
        clicks <- left_join(
            detNameHolder,
            clicks,
            by=c('UID', 'eventId', 'Channel'),
            relationship='many-to-one'
        )
        x <- detDataToStudy(x, clicks)
        return(x)
    }
    clicks
}

# marks keepClick FALSE if moving too fast
# PAMpal util for above prob rename to swim speed?
filterBySwimSpeed <- function(data, time=30, depth=NULL, speed=NULL) {
    if(is.null(depth) && is.null(speed)) {
        stop('Either depth or swim speed (m/s) must be provided')
    }
    if(is.null(depth)) {
        depth <- time * speed
    }
    # init values - keep only those non-NA
    if(!'keepClick' %in% colnames(data)) {
        data$keepClick <- !is.na(data$maxDepth)
    }
    if(all(is.na(data$maxDepth))) {
        return(data)
    }
    if(nrow(data) <= 1) {
        return(data)
    }
    for(i in seq_len(nrow(data))) {
        tDiffs <- abs(as.numeric(difftime(data$UTC[i], data$UTC[-i], units='secs')))
        dDiffs <- abs(data$maxDepth[i] - data$maxDepth[-i])
        # only keep if we wer alrady keeping and this test is still good
        data$keepClick[i] <- data$keepClick[i] &
            isTRUE(any((tDiffs <= time) & (dDiffs <= depth)))
    }
    data
}
