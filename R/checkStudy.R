#' @title Check an AcousticStudy Object for Issues
#'
#' @description Checks for any possible issues in an \linkS4class{AcousticStudy}
#'   object, printing statements to the console
#'
#' @param x an \linkS4class{AcousticStudy} object
#'
#' @return invisibly returns \code{x} after printing possible issues
#'
#' @examples
#'
#' data(exStudy)
#'
#' # currently only checks if any peak frequencies are 0, so well force this
#' exStudy[[1]][[1]]$peak <- 0
#' checkStudy(exStudy)
#'
#' @export
#'
#' @author Taiki Sakai \email{taiki.sakai@@noaa.gov}
#'
checkStudy <- function(x) {
    tryCatch({
        det <- getDetectorData(x)
        if(nrow(det$click) > 0 &&
           'peak' %in% colnames(det$click)) {
            if(any(det$click$peak == 0, na.rm=TRUE)) {
                peak0Msg <- paste0('Some clicks had a peak frequency of 0 Hz,',
                                   ' consider adjusting the filter parameter',
                                   ' or adding a calibration function.')
                warning(peak0Msg, call. = FALSE)
            }
        }
    },
    error = function(e) {
        return(invisible(x))
    })
    invisible(x)
}

# threshold max distance between consec detections 2 hours
#
