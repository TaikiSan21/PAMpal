#' @title Get Binary Data for Detections from an Acoustic Event
#'
#' @description Fetches matching binary data from a single or multiple
#'   detections in an \linkS4class{AcousticEvent} object
#'
#' @param x a \linkS4class{AcousticStudy} object, a list of \linkS4class{AcousticEvent}
#'   objects, or a single \linkS4class{AcousticEvent} object
#' @param UID the UID(s) of the individual detections to fetch the binary
#'   data for
#' @param quiet logical flag to quiet some warnings, used internally
#' @param \dots additional arguments to pass to
#'   \code{\link[PamBinaries]{loadPamguardBinaryFile}}
#'
#' @return a list of \code{PamBinary} objects for each \code{UID}
#'
#' @author Taiki Sakai \email{taiki.sakai@@noaa.gov}
#'
#' @importFrom dplyr bind_rows
#' @importFrom PamBinaries loadPamguardBinaryFile
#' @export
#'
getBinaryData <- function(x, UID, quiet=FALSE, ...) {
    if(is.AcousticStudy(x)) {
        # only look in events that have our UIDs
        # hasUID <- whereUID(x, UID)
        # hasUID <- unique(unlist(hasUID))
        # hasUID <- hasUID[!is.na(hasUID)]
        # x <-events(x)[hasUID]
        x <- events(x)
    }
    if(is.list(x)) {
        # only do AcEv in list, run quietly bc we know not all UIDs in all events
        isAcev <- sapply(x, is.AcousticEvent)
        if(all(!isAcev)) {
            warning('No AcousticEvents found, check inputs.')
            return(NULL)
        }
        if(any(!isAcev)) {
            warning('Not all inputs were AcousticEvent objects.')
        }
        result <- lapply(x[isAcev], function(y) {
            getBinaryData(y, UID, quiet=TRUE, ...)
        })
        names(result) <- NULL
        return(unlist(result, recursive=FALSE))
    }
    if(!is.AcousticEvent(x)) {
        warning('This is not an AcousticEvent object.')
        return(NULL)
    }
    # from here we know its an AcEv
    allBinaries <- files(x)$binaries
    # find matching UID from dets
    bins <- bind_rows(
        lapply(detectors(x), function(df) {
            df[df[['UID']] %in% UID, c('UTC', 'UID', 'BinaryFile')]
        }))
    if(is.null(bins) ||
       nrow(bins) == 0) {
        if(!quiet) {
            warning('No matches found for UID(s) ', paste(UID, collapse=', '), '.')
        }
        return(NULL)
    }
    bins <- unique(bins)

    if(length(settings(x)$sr) == 1) {
        bins$sr <- settings(x)$sr
    } else if(length(settings(x)$sr) > 1) {
        trySr <- matchSR(bins, files(x)$db, safe=TRUE)
        if(is.null(trySr)) {
            warning('Multiple sample rates present for event ', id(x),', but not able to read',
                    ' from database ', files(x)$db,'.\nSample rate will not be attached,',
                    ' check that file exists to fix.')
        } else {
            bins <- trySr
        }
    }

    nIn <- sapply(UID, function(y) {
        sum(y == bins[['UID']])
    })
    if(any(nIn == 0) &&
       !quiet) {
        warning('No matches found for UID(s) ',
                paste(UID[nIn == 0], collapse=', '), '.')
    }
    if(any(nIn > 1)) {
        warning('Multiple matches found for UID(s) ',
                paste(UID[nIn > 1], collapse=', '), '.')
        print(bins[bins[['UID']] %in% UID[nIn > 1],])
    }
    # Bins is detector data
    result <- lapply(unique(bins$BinaryFile), function(bin) {
        # this has full path name
        fullBin <- grep(bin, unique(allBinaries), value = TRUE)
        if(length(fullBin)==0) {
            warning('Binary file ', bin, ' not found in files slot.')
            return(NULL)
        }
        if(length(fullBin) > 1) {
            warning('More than one binary file found for ', bin, ' only using first one.')
            fullBin <- fullBin[1]
        }
        if(!file.exists(fullBin)) {
            warning('Binary file ', fullBin, ' does not exist, please check (this file name',
                    ' should be the full file path).')
            return(NULL)
        }
        data <- loadPamguardBinaryFile(fullBin, skipLarge = FALSE, convertDate = TRUE,
                               keepUIDs = bins[['UID']][bins$BinaryFile == bin], ...)$data
        # just for useful stuff later in other functions and if theres two you see what we chose
        ##### I DONT KNOW IF THIS IS A GOOD IDEAAAAA
        if('sr' %in% names(bins)) {
            for(i in names(data)) {
                data[[i]]$sr <- bins$sr[bins$UID == i][1]
            }
        }
        data
    })
    # list named by UID as result
    unlist(result, recursive = FALSE)
}
