#' @title Constructor for PAMpalSettings Object
#'
#' @description Create a PAMpalSettings object. Any values that are not supplied
#'   will be asked for interactively.
#'
#' @param db the full path to a PamGuard database file
#' @param binaries a list with items "folder" containing the directory of the
#'   PamGuard binary files, and "list" containing the full path to each individual
#'   binary file.
#' @param calibration a named list of calibration functions to apply while
#'   applying functions from the "functions" slot. Should named by the
#'   PamGuard module, same as the "functions"
#'
#' @return A PAMpalSettings object
#'
#' @author Taiki Sakai \email{taiki.sakai@@noaa.gov}
#'
#' @importFrom methods new
#' @export
#'
PAMpalSettings <- function(db=NULL, binaries=NULL, calibration=NULL) {
    pps <- new('PAMpalSettings')
    pps <- addDatabase(pps, db)
    pps <- addBinaries(pps, binaries)
    cat('Default included function(s) are "standardClickCalcs" for the "ClickDetector" module,',
        '"roccaWhistleCalcs" for the "WhistlesMoans" module,',
        'and "standardCepstrumCalcs" for the "Cepstrum" module.\n')
    pps <- addFunction(pps, standardClickCalcs, 'ClickDetector')
    pps <- addFunction(pps, roccaWhistleCalcs, 'WhistlesMoans')
    pps <- addFunction(pps, standardCepstrumCalcs, 'Cepstrum')
    # pps <- addCalibration(pps, calibration)
    pps
}
