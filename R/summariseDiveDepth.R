#' @title Summarise Dive Depth
#'
#' @description Summarise results of dive depth estimation using
#'   \link{calculateEchoDepth} and related functions
#'
#' @param x an \linkS4class{AcousticStudy} that has been
#'   processed with \link{calculateEchoDepth}
#' @param hpDepthError hydrophone depth error to use for error estimation
#' @param locType name of localization, note that this function is not computing
#'   any localization, only using previously calculated
#'
#' @author Taiki Sakai \email{taiki.sakai@@noaa.gov}
#'
#' @return a dataframe with columns summarising the estimated dive depth
#'   for each event in \code{x}
#'
#' @examples
#' # example not run because \link{calculateEchoDepth} must be run first,
#' # and it requires a large amount of data not stored in the package
#' \dontrun{
#' study <- calculateEchoDepth(study, wav='path/to/wavFiles')
#' summariseDiveDepth(study)
#' }
#'
#' @export
#'
summariseDiveDepth <- function(x, hpDepthError=1, locType='PGTargetMotion') {
    clicks <- getClickData(x, measures=FALSE)
    clicks <- distinct(dropCols(clicks, 'detectorName'))
    if(!'maxDepth' %in% colnames(clicks)) {
        stop('Data has not been processed with "calculateEchoDepth"')
    }
    if(!'keepClick' %in% colnames(clicks)) {
        warning('It is recommended to filter with "runDepthReview" first.')
        clicks$keepClick <- !is.na(clicks$maxDepth)
    }
    locData <- bind_rows(lapply(events(x), function(x) {
        localizations(x)[[locType]]
    }), .id='eventId')
    hasLoc <- locData$eventId[!is.na(locData$locLat)]
    clicks <- filter(clicks,
                     .data$eventId %in% hasLoc,
                     .data$keepClick)
    clicks$hpDepthPct <- hpDepthError / clicks$hpDepth
    clicks <- left_join(clicks, locData[c('eventId', 'perpDist', 'perpDistErr')], by='eventId') %>%
        mutate(srPct = .data$perpDistErr / .data$perpDist,
               errPct = sqrt(.data$hpDepthPct^2 + .data$srPct^2),
               depthErr = .data$maxDepth * .data$errPct)

    smry <- clicks %>%
        group_by(.data$eventId) %>%
        summarise(meanDepth = mean(.data$maxDepth),
                  minDepth = min(.data$maxDepth),
                  maxDepth = max(.data$maxDepth),
                  meanErrPct = mean(.data$errPct),
                  originalDistance = mean(.data$perpDist),
                  originalDistErr = mean(.data$perpDistErr)) %>%
        mutate(correctedDistance = suppressWarnings(sqrt(.data$originalDistance^2 - .data$meanDepth^2)))
    smry$correctedDistance[is.na(smry$correctedDistance)] <- 0
    smry
}
