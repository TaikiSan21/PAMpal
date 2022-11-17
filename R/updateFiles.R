#' @title Update Location of Files in an AcousticStudy
#'
#' @description Updates the stored locations of binary, database, and/or recording
#'   files in the \code{files} slots of an \linkS4class{AcousticStudy} and all
#'   \linkS4class{AcousticEvent} objects within. Runs interactively to prompt
#'   users to select folders if missing files are found. Typically used after changing
#'   computers, or if original data was on an external hard drive. If any
#'   missing files are not able to be located, they will be kept in the files
#'   slot so that this function can be run again
#'
#' @param x an \linkS4class{AcousticStudy} or \linkS4class{AcousticEvent}
#'   object
#' @param bin folder containing updated binary file locations. If
#'   \code{NULL} (default), user will be prompted to select a folder.
#'   If \code{NA}, binary files will be skipped.
#' @param db single file or folder containing updated database file locations.
#'   \code{NULL} (default), user will be prompted to select a folder.
#'   If \code{NA}, database files will be skipped.
#' @param recording folder containing updated recording file locations. If
#'   \code{NULL} (default), user will be prompted to select a folder.
#'   If \code{NA}, recording files will be skipped.
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
    # do databases
    # this check is for weirdness in older v of pampal
    if(is.list(files(x)$db)) {
        files(x)$db <- unlist(files(x)$db)
    }
    dbExists <- file.exists(files(x)$db)
    if(all(dbExists)) {
        db <- character(0)
    } else {
        db <- fileLister(db, label = 'database', pattern='sqlite', verbose=verbose)
        updatedDbs <- fileMatcher(files(x)$db, db)
        if(verbose) {
            cat(paste0('Updated the locations of ',
                       sum(!dbExists)-sum(!file.exists(updatedDbs)),
                       ' out of ', sum(!dbExists),
                       ' missing database files.'), '\n')
        }
        files(x)$db <- updatedDbs
    }
    # do binaries
    # this check is for weirdness in older v of pampal
    if(is.list(files(x)$binaries)) {
        files(x)$binaries <- unlist(files(x)$binaries)
    }
    binExists <- file.exists(files(x)$binaries)
    if(all(binExists)) {
        bin <- character(0)
    } else {
        bin <- fileLister(bin, label = 'binary', pattern = ppVars()$binPattern, verbose=verbose)
        updatedBins <- fileMatcher(files(x)$binaries, bin)
        if(verbose) {
            cat(paste0('Updated the locations of ',
                       sum(!binExists) - sum(!file.exists(updatedBins)),
                       ' out of ', sum(!binExists),
                       ' missing binary files.'), '\n')
        }
        files(x)$binaries <- updatedBins
    }
    # do recordings
    if(!is.null(files(x)$recordings$file) &&
       any(!file.exists(files(x)$recordings$file))) {
        fileExists <- file.exists(files(x)$recordings$file)
        recording <- fileLister(recording, label='recording', pattern='wav$', verbose=verbose)
        updatedRecs <- fileMatcher(files(x)$recordings$file, recording)
        if(verbose) {
            cat(paste0('Updated the locations of ',
                       sum(!fileExists) - sum(!file.exists(updatedRecs)),
                       ' out of ', sum(!fileExists),
                       ' missing recording files.'), '\n')

        }
        files(x)$recordings$file <- updatedRecs
    }
    # from here is adjusting events, they dont carry recording info so we can stop here
    # if thats all we changed
    if((length(bin) == 0) &&
       (length(db) == 0)) {
        return(x)
    }
    if(is.AcousticStudy(x)) {
        for(e in seq_along(events(x))) {
            # be quiet for every event, study level should give best summary
            events(x)[[e]] <- updateFiles(events(x)[[e]], bin=bin, db=db, recording=recording, verbose=FALSE)
        }
    }
    x <- .addPamWarning(x)
    x
}

# match and replace
# old is vector or df, new is vector of new file names
fileMatcher <- function(old, new) {
    if(length(new) == 0) {
        return(old)
    }
    if(is.data.frame(old)) {
        old$file <- fileMatcher(old$file, new)
        return(old)
    }
    new <- normalizePath(new)
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

fileLister <- function(x, label, pattern, verbose=TRUE) {
    if(is.null(x)) {
        downLab <- tolower(label)
        upLab <- paste0(toupper(substr(downLab, 1, 1)),
                        substr(downLab, 2, 100))
        cat('Found missing ', downLab, ' files, please select a new folder',
            'where they can be found.\n', sep='')
        # if(!isAvailable('1.1.287')) {
        #     bin <- choose.dir(caption = 'Choose Binary Folder:')
        # } else {
        #     bin <- selectDirectory(caption = 'Choose Binary Folder:', path = getwd())
        # }
        x <- tk_choose.dir(caption = paste0('Choose ', upLab, ' Folder:'), default = getwd())
    }
    x <- x[!is.na(x)]
    if(length(x) == 0) {
        return(character(0))
    }
    if(length(x) == 1 &&
       dir.exists(x)) {
        files <- list.files(x, full.names=TRUE, recursive=TRUE, pattern=pattern)
        if(verbose &&
           length(files) == 0) {
            cat('No ', downLab, ' files found in this folder.\n', sep='')
        }
        return(files)
    }
    DNE <- !file.exists(x)
    if(!any(DNE)) {
        return(x)
    }
    pamWarning('Files ', paste0(x[DNE], collapse = ', '), ' could not be located.')
    x[!DNE]
}
