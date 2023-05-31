#' @title Subsample Detectors in AcousticStudy
#'
#' @description samples either a fraction or fixed number from each detector
#'   in each event of an AcousticStudy
#'
#' @param x \linkS4class{AcousticStudy} object
#' @param n if less than 1, proportion of rows to sample from each detector.
#'   If 1 or more, the number of rows to sample from each detector. If \code{n}
#'   is negative, this proportion or number will be removed from each detector.
#'
#' @details Uses \link[dplyr]{slice_sample} to do the sampling, \code{n}
#'   is converted either to \code{prop} or \code{n} based on its size.
#'   Negative values are treated the same as in \link[dplyr]{slice_sample}
#'
#' @return subsampled AcousticStudy \code{x}
#'
#' @author Taiki Sakai \email{taiki.sakai@@noaa.gov}
#'
#' @examples
#'
#' data(exStudy)
#' nDetections(exStudy)
#' halfData <- sampleDetector(exStudy, n=0.5)
#' # there are odd numbers of rows in some detectors, so less than half
#' nDetections(exStudy)
#' oneDetPerDetector <- sampleDetector(exStudy, n=1)
#' nDetections(exStudy)
#'
#' @importFrom dplyr slice_sample
#' @export
#'
sampleDetector <- function(x, n=1) {
    if(!is.AcousticStudy(x)) {
        stop('"x" must be an AcousticStudy')
    }
    if(abs(n) < 1) {
        sampleFun <- function(d) {
            slice_sample(d, prop=n)
        }
    } else {
        sampleFun <- function(d) {
            slice_sample(d, n=n)
        }
    }
    events(x) <- lapply(events(x), function(e) {
        for(d in seq_along(detectors(e))) {
            ct <- attr(e[[d]], 'calltype')
            e[[d]] <- sampleFun(e[[d]])
            attr(e[[d]], 'calltype') <- ct
        }
        nDets <- sapply(detectors(e), nrow)
        detectors(e) <- detectors(e)[nDets > 0]
        e
    })
    x
}
