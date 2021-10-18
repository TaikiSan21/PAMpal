#' @title Get Raw Binary Data for Detections
#'
#' @description Fetches matching binary data from a single or multiple
#'   detections in an \linkS4class{AcousticEvent} object
#'
#' @param x a \linkS4class{AcousticStudy} object, a list of \linkS4class{AcousticEvent}
#'   objects, or a single \linkS4class{AcousticEvent} object
#' @param UID the UID(s) of the individual detections to fetch the binary
#'   data for
#' @param type detection type 
#' @param quiet logical flag to quiet some warnings, used internally and should generally
#'   not be changed from default \code{FALSE}
#' @param \dots additional arguments to pass to
#'   \code{\link[PamBinaries]{loadPamguardBinaryFile}}
#'
#' @return a list of \code{PamBinary} objects for each \code{UID}
#'
#' @author Taiki Sakai \email{taiki.sakai@@noaa.gov}
#'
#' @examples
#'
#' data(exStudy)
#' binData <- getBinaryData(exStudy, UID = 8000003)
#' # works with multiple UIDs, if UIDs arent present they will be ignored
#' binData <- getBinaryData(exStudy, UID = c(8000003, 529000024, 1))
#'
#' @importFrom dplyr bind_rows
#' @importFrom PamBinaries loadPamguardBinaryFile
#' @export
#'
getBinaryData <- function(x, UID, type=c('click', 'whistle', 'cepstrum'), quiet=FALSE, ...) {
    # if(is.AcousticStudy(x)) {
    #     x <- events(x)
    # }
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
            getBinaryData(y, UID, type=type, quiet=TRUE, ...)
        })
        names(result) <- NULL
        result <- unlist(result, recursive = FALSE)
        # if any repeated UIDs in separate events only keep the first one
        if(length(result) > length(unique(names(result)))) {
            drop <- NULL
            for(i in unique(names(result))) {
                whichId <- which(names(result) == i)
                drop <- c(drop, whichId[-1])
            }
            result <- result[-drop]
        }
        return(result)
    }
    if(!is.AcousticEvent(x) &&
       !is.AcousticStudy(x)) {
        warning('This is not an AcousticEvent or AcousticStudy object.')
        return(NULL)
    }
    # fr
    
    allBinaries <- unique(files(x)$binaries)
    # find matching UID from dets
    bins <- bind_rows(
        lapply(getDetectorData(x), function(df) {
            df[df[['UID']] %in% UID, c('UTC', 'UID', 'BinaryFile', 'detectorName')]
        }))
    # bins$detectorName <- gsub('_[0-9]{0,3}$', '', bins$detectorName)
    typeMatch <- vector('character', length=length(type))
    for(t in seq_along(type)) {
        typeMatch[t] <- switch(type[t],
                            'click' = '^Click_Detector_|^SoundTrap_Click_Detector_',
                            'whistle' = '^WhistlesMoans_',
                            'cepstrum' = '^WhistlesMoans_'
        )
    }
    typeMatch <- paste0(typeMatch, collapse='|')
    bins <- bins[grepl(typeMatch, bins$BinaryFile), ]
    if(is.null(bins) ||
       nrow(bins) == 0) {
        if(!quiet) {
            warning('UID(s) ', printN(UID, 6), ' were not found in data.')
        }
        return(NULL)
    }
    bins <- unique(bins)
    # just doing this bec i goofed earlier and some people had data where this
    # never happened in processing. BinaryFile should already be basename
    bins$BinaryFile <- basename(bins$BinaryFile)
    if(!is.null(getSr(x, type=type, name=bins$detectorName[1], bins$UTC))) {
        bins$sr <- getSr(x, type, bins$detectorName, bins$UTC)
    # } else if(length(settings(x)$sr) == 1) {
    #     bins$sr <- settings(x)$sr
    # } else if(length(settings(x)$sr) > 1) {
    #     trySr <- matchSR(bins, files(x)$db, safe=TRUE)
    #     if(is.null(trySr)) {
    #         warning('Multiple sample rates present for event ', id(x),', but not able to read',
    #                 ' from database ', files(x)$db,'.\nSample rate will not be attached,',
    #                 ' check that file exists to fix.')
    #     } else {
    #         bins <- trySr
    #     }
    # }
    }
    bins$detectorName <- NULL
    bins <- unique(bins)
    nIn <- sapply(UID, function(y) {
        sum(y == bins[['UID']])
    })
    if(any(nIn == 0) &&
       !quiet) {
        warning('No matches found for UID(s) ',
                printN(UID[nIn == 0], 6), '.')
    }
    if(any(nIn > 1)) {
        warning('Multiple matches found for UID(s) ',
                printN(UID[nIn > 1], 6), '.')
        # print(bins[bins[['UID']] %in% UID[nIn > 1],])
    }
    # Bins is detector data
    result <- lapply(unique(bins$BinaryFile), function(bin) {
        # this has full path name
        fullBin <- grep(bin, allBinaries, value = TRUE, fixed=TRUE)
        if(length(fullBin)==0) {
            warning('Binary file ', bin, ' not found in files slot.', call.=FALSE)
            return(NULL)
        }
        data <- list()
        matchUID <- rep(FALSE, length(bins[['UID']][bins$BinaryFile == bin]))
        names(matchUID) <- bins[['UID']][bins$BinaryFile == bin]
        for(i in seq_along(fullBin)) {

            # if(length(fullBin) > 1) {
            #     warning('More than one binary file found for ', bin, ' only using first one.')
            #     fullBin <- fullBin[1]
            # }
            if(!file.exists(fullBin[i])) {
                warning('Binary file ', fullBin[i], ' could not be located, if files have moved ',
                        'please use "updateFiles" function.', call.=FALSE)
                # return(NULL)
                next
            }
            thisData <- loadPamguardBinaryFile(fullBin[i], skipLarge = FALSE, convertDate = TRUE,
                                             keepUIDs = names(matchUID)[!matchUID], ...)$data
            # just for useful stuff later in other functions and if theres two you see what we chose
            ##### I DONT KNOW IF THIS IS A GOOD IDEAAAAA
            if(length(thisData) > 0) {
                matchUID[names(thisData)] <- TRUE
                data <- c(data, thisData)
            }
            if(all(matchUID)) {
                break
            }
        }
        if(length(data) == 0) {
            warning('No data found for binary file ', bin, call.=FALSE)
            return(NULL)
        }
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
