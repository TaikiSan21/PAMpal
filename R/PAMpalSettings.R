#' @title Constructor for PAMpalSettings Object
#'
#' @description Create a PAMpalSettings object. Any values that are not supplied
#'   will be asked for interactively. Three processing functions will also be added
#'   by default: \link{standardClickCalcs}, \link{roccaWhistleCalcs}, and
#'   \link{standardCepstrumCalcs}
#'
#' @param db the full path to a Pamguard database file or folder of databases
#' @param binaries a folder containing Pamguard binary files, all subfolders will
#'   also be added
#' @param settings an XML settings file from Pamguard
#' @param functions a named list of additional functions to add
#' @param verbose logical flag to show messages
#' @param default logical flag to use default measurement function parameters
#' @param \dots values to pass on to default \link{standardClickCalcs} function
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
PAMpalSettings <- function(db=NULL, binaries=NULL, settings=NULL, functions=NULL, verbose=TRUE, default=FALSE, ...) {
    pps <- new('PAMpalSettings')
    pps <- addDatabase(pps, db, verbose)
    pps <- addBinaries(pps, binaries, verbose)
    if(verbose) {
        cat('Default included function(s) are "standardClickCalcs" for the "ClickDetector" module,',
            '"roccaWhistleCalcs" for the "WhistlesMoans" and "GPLDetector" modules,',
            'and "standardCepstrumCalcs" for the "Cepstrum" module.\n')
    }
    pps <- addFunction(pps, standardClickCalcs, 'ClickDetector', verbose=verbose, default=default, ...)
    pps <- addFunction(pps, roccaWhistleCalcs, 'WhistlesMoans', verbose=verbose)
    pps <- addFunction(pps, standardCepstrumCalcs, 'Cepstrum', verbose=verbose)
    pps <- addFunction(pps, roccaWhistleCalcs, 'GPLDetector', verbose=verbose)
    if(!is.null(settings)) {
        pps <- addSettings(pps, settings, verbose=verbose)
    }
    if(!is.null(functions)) {
        fchar <- deparse(substitute(functions))
        fchar <- gsub('^list\\(|\\)$]', '', fchar)
        fchar <- strsplit(fchar, ',')[[1]]
        for(f in seq_along(functions)) {
            thisFchar <- fchar[f]
            thisFchar <- strsplit(thisFchar, '=')[[1]][2]
            thisFchar <- gsub('\\s|\\)$', '', thisFchar)
            attr(functions[[f]], 'fname') <- thisFchar
            pps <- addFunction(pps, functions[[f]], names(functions)[f], verbose=verbose, default=default, ...)
        }
    }
    pps
}
