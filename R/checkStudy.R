#' @title Check an AcousticStudy Object for Issues
#'
#' @description Checks for any possible issues in an \linkS4class{AcousticStudy}
#'   object, printing statements to the console
#'
#' @param x an \linkS4class{AcousticStudy} object
#'
#' @return returns \code{x} after printing possible issues
#'
#' @export
#'
#' @author Taiki Sakai \email{taiki.sakai@@noaa.gov}
#'
checkStudy <- function(x) {
    det <- getDetectorData(x)
    if(nrow(det$click) > 0 &&
       'peak' %in% colnames(det$click)) {
        if(any(det$click$peak == 0)) {
            peak0Msg <- paste0('Some clicks had a peak frequency of 0 Hz,',
                               ' consider adjusting the filter parameter',
                               ' or adding a calibration function.')
            cat('\n', peak0Msg, sep='')
        }
    }
    x
}
