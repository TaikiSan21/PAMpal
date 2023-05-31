#' @title Subsample Detectors in AcousticStudy
#'
#' @description samples either a fraction or fixed number from each detector
#'   in each event of an AcousticStudy
#'
#' @param x \linkS4class{AcousticStudy} object
#' @param n if less than 1, proportion of rows to sample from each detector.
#'   If 1 or more, the number of rows to sample from each detector.
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
    if(is.AcousticEvent(x)) {
        return(sampleEvent(x, n=n))
    }
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
        sampleEvent(e, FUN=sampleFun)
    })
    x
}

sampleEvent <- function(x, FUN=NULL, n=1) {
    if(!is.AcousticEvent(x)) {
        warning('"x" must be an AcousticEvent')
        return(x)
    }
    if(is.null(FUN)) {
        if(abs(n) < 1) {
            FUN <- function(d) {
                slice_sample(d, prop=n)
            }
        } else {
            FUN <- function(d) {
                slice_sample(d, n=n)
            }
        }
    }
    for(d in seq_along(detectors(x))) {
        ct <- attr(x[[d]], 'calltype')
        x[[d]] <- FUN(x[[d]])
        attr(x[[d]], 'calltype') <- ct
    }
    nDets <- sapply(detectors(x), nrow)
    detectors(x) <- detectors(x)[nDets > 0]
    x
}
