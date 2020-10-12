#' @title Calculate a Set of Measurements from a Cepstrum Contour
#'
#' @description Calculate a set of measurements from a cepstrum contour. This
#'   is currently used to measure the inter-click interval of the burst pulse
#'   type calls
#'
#' @param data a list that must have \code{quefrency} the "quefrency" at each
#'   cepstrum contour, \code{sr} the sample rate of the data,
#'   and \code{time} the time in seconds at each bin
#'
#' @return A list with inter-click interval (ici), duration (seconds), and
#'   slope of the ici
#'
#' @author Taiki Sakai \email{taiki.sakai@@noaa.gov}
#'
#' @examples
#'
#' data(testCeps)
#' standardCepstrumCalcs(testCeps)
#'
#' @importFrom stats lm
#' @export
#'
standardCepstrumCalcs <- function(data) {
    neededVals <- c('quefrency', 'time', 'sr')
    missingVals <- neededVals[!(neededVals %in% names(data))]
    if(length(missingVals) > 0) {
        warning('Values for', paste(missingVals, collapse=', '), 'are missing.',
             'These are required for Cepstrum Calculations, please fix.')
        return(NULL)
    }
    nSlices <- length(data$quefrency)
    result <- list(ici = median(data$quefrency/data$sr))
    result$duration <- max(data$time) - min(data$time)
    result$iciSlope <- unname(lm(quefrency ~ time, data=data[c('quefrency', 'time')])$coefficients[2] / data$sr)
    result
}
