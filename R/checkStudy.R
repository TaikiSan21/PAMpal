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
#' # setting up example data
#' exPps <- new('PAMpalSettings')
#' exPps <- addDatabase(exPps, system.file('extdata', 'Example.sqlite3', package='PAMpal'))
#' exPps <- addBinaries(exPps, system.file('extdata', 'Binaries', package='PAMpal'))
#' exClick <- function(data) {
#'     standardClickCalcs(data, calibration=NULL, filterfrom_khz = 0)
#' }
#' exPps <- addFunction(exPps, exClick, module = 'ClickDetector')
#' exPps <- addFunction(exPps, roccaWhistleCalcs, module='WhistlesMoans')
#' exPps <- addFunction(exPps, standardCepstrumCalcs, module = 'Cepstrum')
#' exData <- processPgDetections(exPps, mode='db')
#'
#' # currently only checks if any peak frequencies are 0, so well force this
#' exData[[1]][[1]]$peak <- 0
#' checkStudy(exData)
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
            if(any(det$click$peak == 0)) {
                peak0Msg <- paste0('Some clicks had a peak frequency of 0 Hz,',
                                   ' consider adjusting the filter parameter',
                                   ' or adding a calibration function.')
                warning(peak0Msg, call. = FALSE)
            }
        }
    },
    error = function(e) {
        cat('Please show Taiki this message:\n')
        print(det$click$peak)
        return(invisible(x))
    })
    invisible(x)
}
