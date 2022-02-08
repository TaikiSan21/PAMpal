#' @title Add a Database to a PAMpalSettings Object
#'
#' @description Adds a new function to the "function" slot in a PAMpalSettings
#'   object. Interactively asks for database files if none are supplied as input
#'
#' @param pps a \linkS4class{PAMpalSettings} object to add a database to
#' @param db a database to add
#' @param verbose logical flag to show messages
#'
#' @return the same \linkS4class{PAMpalSettings} object as pps, with the database
#'   \code{db} added to the "db" slot
#'
#' @examples
#'
#' # not recommended to create a pps like this, for example only
#' pps <- new('PAMpalSettings')
#' db <- system.file('extdata', 'Example.sqlite3', package='PAMpal')
#' pps <- addDatabase(pps, db)
#' pps
#'
#' @author Taiki Sakai \email{taiki.sakai@@noaa.gov}
#'
#' @importFrom tcltk tk_choose.files
#' @export
#'
addDatabase <- function(pps, db=NULL, verbose=TRUE) {
    if(is.PAMpalSettings(db)) {
        db <- db@db
    }
    if(is.null(db)) {
        cat('Please select a database file,',
            ' multiple selections are ok..\n')
        db <- tk_choose.files(caption='Select database(s):')
    }
    # Case when cancelled or some weirdness
    if(length(db) == 0) return(pps)

    exists <- file.exists(db)
    if(any(!exists)) {
        warning(paste0('Database(s) ',
                       paste0(db[!exists], collapse=', ')
                       , ' do(es) not exist.'))
        db <- db[exists]
    }
    isSqlite <- grepl('\\.sqlite3$', db)
    if(any(!isSqlite)) {
        warning('Some files selected that are not sqlite3 databases,',
                ' these files have been removed from the selection: ',
                paste0(db[!isSqlite], collapse = ', '))
        db <- db[isSqlite]
    }
    db <- normalizePath(db)
    if(verbose) {
        cat(paste0('Adding ', length(db), ' databases:\n  ', printN(basename(db), 6, collapse='\n  '), '\n'))
        # if(length(db) > 6) {
        #     dbMsg <- paste0(c(basename(db[1:6]),paste0('... (', length(db)-6, ' more not shown)')), collapse = '\n  ')
        # } else {
        #     dbMsg <- paste0(basename(db), collapse=', ')
        # }
        # cat
    }
    pps@db <- unique(c(pps@db, db))
    pps
}
