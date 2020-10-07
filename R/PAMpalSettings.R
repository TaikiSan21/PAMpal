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
#' @param verbose logical flag to show messages
#'
#' @return A PAMpalSettings object
#'
#' @author Taiki Sakai \email{taiki.sakai@@noaa.gov}
#'
#' @examples
#'
#' # can be run with no arguments with popup menu selections
#' if(interactive()) pps <- PAMpalSettings()
#' db <- system.file('extdata', 'Example.sqlite3', package='PAMpal')
#' bin <- system.file('extdata', 'Binaries', package='PAMpal')
#' # or data folders can be supplied ahead of time
#' if(interactive()) pps <- PAMpalSettings(db=db, binaries=bin)
#'
#' @importFrom methods new
#' @export
#'
PAMpalSettings <- function(db=NULL, binaries=NULL, calibration=NULL, verbose=TRUE) {
    pps <- new('PAMpalSettings')
    pps <- addDatabase(pps, db, verbose)
    pps <- addBinaries(pps, binaries, verbose)
    if(verbose) {
        cat('Default included function(s) are "standardClickCalcs" for the "ClickDetector" module,',
            '"roccaWhistleCalcs" for the "WhistlesMoans" module,',
            'and "standardCepstrumCalcs" for the "Cepstrum" module.\n')
    }
    pps <- addFunction(pps, standardClickCalcs, 'ClickDetector', verbose=verbose)
    pps <- addFunction(pps, roccaWhistleCalcs, 'WhistlesMoans', verbose=verbose)
    pps <- addFunction(pps, standardCepstrumCalcs, 'Cepstrum', verbose=verbose)
    # pps <- addCalibration(pps, calibration)
    pps
}
