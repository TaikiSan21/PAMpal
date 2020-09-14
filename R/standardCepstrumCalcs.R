#' @title Calculate a Set of Measurements from a Cepstrum Contour
#'
#' @description Calculate a set of measurements from a cepstrum contour. This
#'   is currently used to measure the inter-pulse interval of the burst pulse
#'   type calls, please refer to JASA PAPER DOOOD for details.
#'
#' @param data a list that must have \code{quefrency} the "quefrency" at each
#'   cepstrum contour, \code{sr} the sample rate of the data,
#'   and \code{time} the time in seconds at each bin
#'
#' @return A list
#'
#' @author Taiki Sakai \email{taiki.sakai@@noaa.gov}
#'
#' @export
#'
standardCepstrumCalcs <- function(data) {
    nSlices <- length(data$quefrency)
    result <- list(ici = median(data$quefrency/data$sr))
    result$duration <- max(data$time) - min(data$time)
    result$avgSlope <- (data$quefrency[nSlices] - data$quefrency[1]) / result$duration
    result
}
