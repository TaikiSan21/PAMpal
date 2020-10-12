#' @title Constructor for PAMpalSettings Object
#'
#' @description Create a PAMpalSettings object. Any values that are not supplied
#'   will be asked for interactively. Three processing functions will also be added
#'   by default: \link{standardClickCalcs}, \link{roccaWhistleCalcs}, and
#'   \link{standardCepstrumCalcs}
#'
#' @param db the full path to a Pamguard database file
#' @param binaries a folder containing Pamguard binary files, all subfolders will
#'   also be added
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
PAMpalSettings <- function(db=NULL, binaries=NULL, verbose=TRUE) {
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
    pps
}
