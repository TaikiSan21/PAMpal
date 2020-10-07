#' @docType data
#' @name exStudy
#' @title Example AcousticStudy Object
#' @description An example AcousticStudy object created using the example
#' PAMpalSettings object provided with the package. Processed with mode='db'
#' @usage data(exStudy)
#' @format a \linkS4class{AcousticStudy} object containing two
#'   \linkS4class{AcousticEvent} objects
#' @keywords datasets
NULL

#' @docType data
#' @name testClick
#' @title A two-channel recording of a delphinid click
#' @description An example delphinid click waveform. This is a two-channel recording
#' of some kind of delphinid click, recorded at 192kHz. There are
#' 800 samples recorded on each channel.
#' @usage data(testClick)
#' @format A list with two items:
#' \describe{
#'   \item{wave}{a matrix with two columns of 800 samples, each column
#'   is a separate recording channel}
#'   \item{sr}{the sample rate of the recording}
#' }
#' @source Southwest Fisheries Science Center / NMFS / NOAA
#' @keywords datasets
NULL

#' @docType data
#' @name testWhistle
#' @title A fake whistle contour
#' @description A manually created fake whistle contour reanging from 1kHz to 3.1kHz
#' @usage data(testWhistle)
#' @format A list with two items:
#' \describe{
#'   \item{freq}{a vector of the frequency contour values in hertz}
#'   \item{time}{a vector of the time values of the contour in seconds}
#' }
#' @keywords datasets
NULL

#' @docType data
#' @name testCeps
#' @title A fake cepstrum contour
#' @description A manually created fake cepstrum contour, mimicing what the output
#'   would be from the Pamguard module and fed into the cepstrum calcs
#' @usage data(testCeps)
#' @format A list with three items:
#' \describe{
#'   \item{quefrency}{a vector of the cepstrum contour bin numbers, not actually quefrency}
#'   \item{time}{a vector of the time values of the cepsturm contour in seconds}
#'   \item{sr}{the sample rate of the recording}
#' }
#' @keywords datasets
NULL
