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
#' @param check logical flag to do extra checking. You do not need to set this
#'   parameter, used internally for speed purposes so certain checks are not
#'   repeated
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
updateFiles <- function(x, bin=NULL, db=NULL, recording=NULL, verbose=TRUE, check=TRUE) {
    # do databases
    # this check is for weirdness in older v of pampal
    if(is.list(files(x)$db)) {
        files(x)$db <- unlist(files(x)$db)
    }
    # if nothing just set this true to avoid error and cause skip in next check
    if(is.null(files(x)$db)) {
        dbExists <- TRUE
    } else {
        dbExists <- file.exists(files(x)$db)
    }
    if(all(dbExists)) {
        db <- character(0)
    } else {
        db <- fileLister(db, label = 'database', pattern='sqlite', verbose=verbose, check=check)
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
    if(is.null(files(x)$binaries)) {
        binExists <- TRUE
    } else {
        binExists <- file.exists(files(x)$binaries)
    }
    if(all(binExists)) {
        bin <- character(0)
    } else {
        bin <- fileLister(bin, label = 'binary', pattern = ppVars()$binPattern, verbose=verbose, check=check)
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
        recording <- fileLister(recording, label='recording', pattern='wav$', verbose=verbose, check=check)
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
        if(verbose) {
            cat('Updating files in events...\n')
            pb <- txtProgressBar(min=0, max=length(events(x)), style=3)
        }
        for(e in seq_along(events(x))) {
            # be quiet for every event, study level should give best summary
            events(x)[[e]] <- updateFiles(events(x)[[e]], bin=bin, db=db, recording=recording, verbose=FALSE, check=FALSE)
            if(verbose) {
                setTxtProgressBar(pb, value=e)
            }
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
    # new <- normalizePath(new, winslash = '/')
    # old <- normalizePath(old, winslash = '/')
    # old <- normalizePath(old, winslash = '/', mustWork = FALSE)
    old <- gsub('\\\\', '/', old)
    new  <- gsub('\\\\', '/', new)
    ###
    toCheck <- rep(TRUE, length(old))
    nCheck <- sum(toCheck)
    checkPos <- 1
    while(nCheck > 0) {
        # this is only incrementing when things dont work
        # so just quit eventually so we dont get stuck
        if(checkPos > 5) {
            break
        }
        pathDiff <- getOnePathDiff(old[toCheck][checkPos], new)
        # case when somehow no match
        if(length(pathDiff) == 1 && is.na(pathDiff)) {
            checkPos <- checkPos + 1
            next
        }
        old[toCheck] <- gsub(pathDiff[1], pathDiff[2], old[toCheck], fixed=TRUE)
        newToCheck <- !old %in% new
        newNCheck <- sum(newToCheck)
        # if we didnt fix any, move one down the line of ones to check
        # and try using that base directory. If we do this too many times
        # then just give up
        if(newNCheck == nCheck) {
            checkPos <- checkPos + 1
        }
        nCheck <- newNCheck
    }
    if(nCheck == 0) {
        return(old)
    }

    # cat(nReps, 'different base directories used')
    # old
    ###
    oldBase <- basename(old)
    newBase <- basename(new)
    newPoss <- new[newBase %in% oldBase]
    newBase <- basename(newPoss)

    matchFun <- function(x) {
        poss <- which(newBase %in% basename(x))
        if(length(poss) == 0) {
            return(NA)
        }
        if(length(poss) == 1) {
            return(poss)
        }
        dirMatch <- checkNextDir(x, newPoss[poss]) # return ix of poss
        poss[dirMatch]
    }

    repIx <- unlist(purrr::map(old, matchFun))
    old[!is.na(repIx)] <- newPoss[repIx[!is.na(repIx)]]
    notIn <- !old %in% new
    nNot <- sum(notIn)
    if(nNot > 0) {
        warning(nNot, ' files could not be updated to new locations')
    }
    old
}

# get path difference for a single old file
getOnePathDiff <- function(old, new) {
    oldBase <- basename(old)
    if(is.null(names(new))) {
        newBase <- basename(new)
    } else {
        newBase <- names(new)
    }
    poss <- which(newBase %in% oldBase)
    if(length(poss) == 0) {
        return(NA)
    }
    if(length(poss) == 1) {
        matchIx <- poss
    }
    if(length(poss) > 1) {
        dirMatch <- checkNextDir(old, new[poss]) # return ix of poss
        matchIx <- poss[dirMatch]
    }
    findPathDiff(old, new[matchIx])
}

# if basenames match, step down through directories until you find actual match
# eg a/b/c/d.csv matches z/x/c/d.csv and z/x/y/d.csv both at first so find best
checkNextDir <- function(x,y) {
    if(x == '.') {
        return(1)
    }
    poss <- which(basename(dirname(y)) %in% basename(dirname(x)))
    if(length(poss) == 0) {
        return(1)
    }
    if(length(poss) == 1) {
        return(poss)
    }
    checkNextDir(dirname(x), dirname(y))
}

# get base parts of paths that differ - c/taiki/files vs c/greg/files
# use this to gsub c/taiki for c/greg
findPathDiff <- function(x, y) {
    # need to check for \\\\ start
    xhas4 <- identical(substr(x, 1, 2), '\\\\')
    yhas4 <- identical(substr(y, 1, 2), '\\\\')
    if(xhas4) {
        x <- substr(x, 3, nchar(x))
    }
    if(yhas4) {
        y <- substr(y, 3, nchar(y))
    }
    splitx <- strsplit(x, '/|\\\\')[[1]]
    splity <- strsplit(y, '/|\\\\')[[1]]
    nx <- length(splitx)
    ny <- length(splity)
    for(i in 1:(max(nx, ny))) {
        if(splitx[nx-i+1] != splity[ny-i+1]) {
            break
        }
    }
    outx <- paste0(splitx[1:(nx-i+1)], collapse='/')
    outy <- paste0(splity[1:(ny-i+1)], collapse='/')
    if(xhas4) {
        outx <- paste0('\\\\', outx)
    }
    if(yhas4) {
        outy <- paste0('\\\\', outy)
    }
    c(outx, outy)
}

fileLister <- function(x, label, pattern, verbose=TRUE, check=TRUE) {
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
        files <- normalizePath(files, winslash = '/')
        names(files) <- basename(files)
        return(files)
    }
    if(is.null(names(x))) {
        names(x) <- basename(x)
    }
    if(isFALSE(check)) {
        return(x)
    }
    DNE <- !file.exists(x)
    if(!any(DNE)) {
        return(x)
    }
    pamWarning('Files ', paste0(x[DNE], collapse = ', '), ' could not be located.')
    x[!DNE]
}
