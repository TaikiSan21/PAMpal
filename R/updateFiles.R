#' @title Update Location of Binary and Databases Files
#'
#' @description Updates the stored locations of binary and database files in
#'   the \code{files} slots of all \linkS4class{AcousticStudy} and
#'   \linkS4class{AcousticEvent} objects. Typically used after changing
#'   computers, or if original data was on an external hard drive. If any
#'   missing files are not able to be located, they will be kept in the files
#'   slot so that this function can be run again
#'
#' @param x an \linkS4class{AcousticStudy} or \linkS4class{AcousticEvent}
#'   object
#' @param bin folder containing updated binary file locations. If
#'   \code{NULL} (default), user will be prompted to select a folder
#' @param db folder containing updated database file locations. If
#'   \code{NULL} (default), user will be prompted to select a folder
#' @param quiet logical flag to print messages about success of replacement
#'
#' @return the same \linkS4class{AcousticStudy} and
#'   \linkS4class{AcousticEvent} object as \code{x} with
#'   updated file locations
#'
#' @author Taiki Sakai \email{taiki.sakai@@noaa.gov}
#'
#' @export
#'
updateFiles <- function(x, bin=NULL, db=NULL, quiet=FALSE) {
    dbExists <- file.exists(files(x)$db)
    if(any(!dbExists)) {
        if(is.null(db)) {
            cat('Found missing databases, please select a new folder',
                'where they can be found.\n')
            # if(!isAvailable('1.1.287')) {
            #     db <- choose.dir(caption = 'Choose Database Folder:')
            # } else {
            #     db <- selectDirectory(caption = 'Choose Database Folder:', path = getwd())
            # }
            db <- tk_choose.dir(caption = 'Choose Database Folder:', default = getwd())
        }
        newDbs <- list.files(db, full.names=TRUE, pattern='sqlite')
        updatedDbs <- fileMatcher(files(x)$db, newDbs)
        if(!quiet) {
            if(length(newDbs) == 0) {
                cat('No database files found in this folder.\n')
            }
            cat(paste0('Updated the locations of ',
                       sum(!dbExists)-sum(!file.exists(updatedDbs)),
                       ' out of ', sum(!dbExists),
                       ' missing database files.'), '\n')
        }
        files(x)$db <- updatedDbs
    }
    binExists <- file.exists(files(x)$binaries)
    if(any(!binExists)) {
        if(is.null(bin)) {
            cat('Found missing binary files, please select a new folder',
                'where they can be found.\n')
            # if(!isAvailable('1.1.287')) {
            #     bin <- choose.dir(caption = 'Choose Binary Folder:')
            # } else {
            #     bin <- selectDirectory(caption = 'Choose Binary Folder:', path = getwd())
            # }
            bin <- tk_choose.dir(caption = 'Choose Binary Folder:', default = getwd())
        }
        newBins <- list.files(bin, recursive = TRUE, full.names = TRUE, pattern ='(Clicks|WhistlesMoans).*pgdf$')

        updatedBins <- fileMatcher(files(x)$binaries, newBins)
        if(!quiet) {
            if(length(newBins) == 0) {
                cat('No binaray files found in this folder.\n')
            }
            cat(paste0('Updated the locations of ',
                       sum(!binExists) - sum(!file.exists(updatedBins)),
                       ' out of ', sum(!binExists),
                       ' missing binary files.'), '\n')
        }
        files(x)$binaries <- updatedBins
    }
    if(is.AcousticStudy(x)) {
        for(e in seq_along(events(x))) {
            events(x)[[e]] <- updateFiles(events(x)[[e]], bin=bin, db=db, quiet=TRUE)
        }
    }
    x
}

# match and replace
fileMatcher <- function(old, new) {
    if(length(new) == 0) {
        return(old)
    }
    hasReplacement <- sapply(old, function(x) any(grepl(basename(x), basename(new))))
    isReplacement <- sapply(new, function(x) any(grepl(basename(x), basename(old))))
    c(old[!hasReplacement], new[isReplacement])
}
