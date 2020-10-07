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
#' @param recording folder containing updated recording file locations. If
#'   \code{NULL} (default), user will be prompted to select a folder
#' @param verbose logical flag to print messages about success of replacement
#'
#' @return the same \linkS4class{AcousticStudy} and
#'   \linkS4class{AcousticEvent} object as \code{x} with
#'   updated file locations
#'
#' @author Taiki Sakai \email{taiki.sakai@@noaa.gov}
#'
#' @examples
#'
#' data(exStudy)
#' # files in exStudy will have paths local to package creator's computer
#' files(exStudy)$db
#' file.exists(files(exStudy)$db)
#' files(exStudy)$binaries
#' file.exists(files(exStudy)$binaries)
#' # folder with example DB
#' db <- system.file('extdata', package='PAMpal')
#' # folder with example binaries
#' bin <- system.file('extdata', 'Binaries', package='PAMpal')
#' exStudy <- updateFiles(exStudy, db=db, bin=bin)
#' files(exStudy)$db
#' file.exists(files(exStudy)$db)
#' files(exStudy)$binaries
#' file.exists(files(exStudy)$binaries)
#'
#' @export
#'
updateFiles <- function(x, bin=NULL, db=NULL, recording=NULL, verbose=TRUE) {
    dbExists <- file.exists(files(x)$db)
    if(all(dbExists)) {
        db <- NA
    }
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
        if(is.na(db)) {
            newDbs <- character(0)
        } else {
            newDbs <- list.files(db, full.names=TRUE, pattern='sqlite')
        }
        updatedDbs <- fileMatcher(files(x)$db, newDbs)
        if(verbose) {
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
    if(all(binExists)) {
        bin <- NA
    }
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
        if(is.na(bin)) {
            newBins <- character(0)
        } else {
            newBins <- list.files(bin, recursive = TRUE, full.names = TRUE, pattern ='(Clicks|WhistlesMoans).*pgdf$')
        }

        updatedBins <- fileMatcher(files(x)$binaries, newBins)
        if(verbose) {
            if(length(newBins) == 0) {
                cat('No binary files found in this folder.\n')
            }
            cat(paste0('Updated the locations of ',
                       sum(!binExists) - sum(!file.exists(updatedBins)),
                       ' out of ', sum(!binExists),
                       ' missing binary files.'), '\n')
        }
        files(x)$binaries <- updatedBins
    }
    if(!is.null(files(x)$recordings$file) &&
       any(!file.exists(files(x)$recordings$file))) {
        fileExists <- file.exists(files(x)$recordings$file)
        if(is.null(recording)) {
            cat('Found missing recording files, please select a new folder',
                'where they can be found.\n')
            recording <- tk_choose.dir(caption='Choose Recording Folder:', default=getwd())
        }
        if(is.na(recording)) {
            newRecs <- character(0)
        } else{
            newRecs <- list.files(recording, recursive=TRUE, full.names=TRUE, pattern='wav$')
        }
        updatedRecs <- fileMatcher(files(x)$recordings$file, newRecs)
        if(verbose) {
            if(length(newRecs) == 0) {
                cat('No recording files found in this folder.\n')
            }
            cat(paste0('Updated the locations of ',
                       sum(!fileExists) - sum(!file.exists(updatedRecs)),
                       ' out of ', sum(!fileExists),
                       ' missing recording files.'), '\n')

        }
        files(x)$recordings$file <- updatedRecs
    }

    if(is.AcousticStudy(x)) {
        for(e in seq_along(events(x))) {
            if(is.na(bin) && is.na(db)) next
            # be quiet for every event, study level should give best summary
            events(x)[[e]] <- updateFiles(events(x)[[e]], bin=bin, db=db, recording=recording, verbose=FALSE)
        }
    }
    x
}

# match and replace
fileMatcher <- function(old, new) {
    if(length(new) == 0) {
        return(old)
    }
    if(is.data.frame(old)) {
        old$file <- fileMatcher(old$file, new)
        return(old)
    }
    repIx <- sapply(old, function(x) {
        rix <- which(grepl(basename(x), basename(new)))
        if(length(rix) == 0) {
            return(0)
        }
        rix[1]
    })
    old[repIx > 0] <- new[repIx[repIx > 0]]
    old
}
