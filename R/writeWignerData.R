#' @title Write Wigner Transform Data of Click Detections to Disk
#'
#' @description Create Wigner-Ville transform data of click clips from all
#'   detections and save them to disk. A CSV file will also be written that
#'   lists all UIDs contained in the output
#'
#' @param x \linkS4class{AcousticStudy} object containing data to make Wigner data for
#' @param n number of frequency bins for Wigner transform (recommended power of 2)
#' @param t number of samples to use for the click clip passed to the transform
#' @param outDir directory to write data to
#' @param mode specifies the kind of output that will be created, currently only
#'   supports creating NumPy arrays using the \code{reticulate} package, in future
#'   will support image creation
#' @param progress logical flag to show progress bar
#' @param \dots optional arguments to pass
#'
#' @return A list with two items: \code{files} - a vector of file names
#'   for the Wigner data that were successfully created, any that were not
#'   able to be written will be \code{NA}, and \code{warnings}, a list with
#'   items containing event IDs that triggered any warnings
#'
#' @author Taiki Sakai \email{taiki.sakai@@noaa.gov}
#'
#' @examples
#'
#' data(exStudy)
#' exStudy <- setSpecies(exStudy, method='pamguard')
#' \dontrun{
#' # not running because files are written to disk
#' wigFiles <- writeWignerData(exStudy, outDir = 'WigFolder')
#' }
#'
#' @importFrom reticulate py_install install_miniconda np_array import py_available py_module_available
#' @importFrom PAMmisc wignerTransform
#' @importFrom utils write.csv
#' @export
#'
writeWignerData <- function(x, n=256, t=300, outDir='.', mode='nparray', progress=TRUE, ...) {
    if(!dir.exists(outDir)) {
        dir.create(outDir)
    }
    if(mode == 'nparray') {
        if(!py_available(TRUE)) {
            message('Python is required but could not be found, is it okay to try and install via miniconda?')
            pChoose <- menu(choices = c('Yes', 'No'), title = 'Install Python?')
            if(pChoose != 1) {
                stop('Cannot continue without Python')
            }
            install_miniconda()
        }
        if(!py_module_available('numpy')) {
            py_install('numpy')
        }
        np <- import('numpy')
        wigNpArray <- function(data, n=256, t=300, dir, filename, intMax = 255) {
            if(!dir.exists(dir)) {
                dir.create(dir)
            }
            result <- array(NA, dim = c(n, t, length(data) * 2))
            isNA <- integer(0)
            for(d in seq_along(data)) {
                if(is.null(data[[d]]$wave)) {
                    isNA <- c(isNA, 2*(d-1) + 1:2)
                    next
                }
                for(c in 1:ncol(data[[d]]$wave)) {
                    thisWav <- clipAroundPeak(data[[d]]$wave[, c], t)
                    wig <- wignerTransform(thisWav, n=n, sr=data[[d]]$sr)$tfr
                    wig <- (wig - min(wig)) / (max(wig) - min(wig)) * intMax
                    result[,,2*(d-1) + c] <- wig
                }
            }
            UID <- rep(names(data), each=2)
            if(length(isNA) > 0) {
                result <- result[,,-isNA]
                UID <- UID[-isNA]
            }
            if(dim(result)[3] == 0) {
                return(list(file=NA,
                            UID=NA))
            }
            # UID <- unique(UID)
            # result <- (result - min(result)) / (max(result)-min(result)) * intMax
            result <- np_array(result, dtype='uint8')
            np$savez_compressed(file.path(dir, filename), result)
            list(file=paste0(file.path(dir, filename), '.npz'),
                 UID = UID)
        }
    }

    naSp <- character(0)
    noDet <- character(0)
    noBin <- character(0)
    allFiles <- character(0)
    spList <- unique(species(x))
    badDir <- sapply(spList, function(x) {
        thisDir <- file.path(outDir, x)
        if(dir.exists(thisDir)) {
            return(FALSE)
        } else {
            !dir.create(thisDir, showWarnings=FALSE)
        }
    })
    if(any(badDir)) {
        warning('Unable to create directories for species names ', printN(spList[badDir], 20),
                ', most likely due to symbols in the names. These will be excluded from the output',
                ', change the species names using "setSpecies" with method="reassign" and run again',
                ' for these species.')
        # keepSpec <- spList[!badDir]
        # this doesnt work because keepSpec doesnt exist in that environment, need to evaluate
        # first somehow? FML
        # x <- filter(x, species %in% keepSpec)
        keepSpec <- species(x) %in% spList[!badDir]
        x <- x[keepSpec]
    }

    if(progress) {
        cat('Writing wigner data...\n')
        pb <- txtProgressBar(min=0, max=length(events(x)), style=3)
    }
    for(e in seq_along(events(x))) {
        thisEv <- x[[e]]
        thisDet <- getDetectorData(thisEv)$click
        if(is.null(thisDet) ||
           nrow(thisDet) == 0)  {
            noDet <- c(noDet, id(thisEv))
            if(progress) {
                setTxtProgressBar(pb, value=e)
            }
            next
        }
        if(is.na(species(thisEv)$id)) {
            naSp <- c(naSp, id(thisEv))
            if(progress) {
                setTxtProgressBar(pb, value=e)
            }
            next
        }
        binData <- getBinaryData(thisEv, UID = thisDet$UID)
        if(is.null(binData) ||
           length(binData) == 0) {
            noBin <- c(noBin, id(thisEv))
            if(progress) {
                setTxtProgressBar(pb, value=e)
            }
            next
        }
        wig <- wigNpArray(binData, n=n, t=t, filename = id(thisEv), dir=file.path(outDir, species(thisEv)$id), ...)
        if(is.na(wig$file)) {
            noDet <- c(noDet, id(thisEv))
            if(progress) {
                setTxtProgressBar(pb, value=e)
            }
            next
        }
        allFiles <- c(allFiles, wig$file)
        df <- data.frame(UID=wig$UID, event=id(thisEv), species = species(thisEv)$id, stringsAsFactors = FALSE)
        write.csv(df, file=paste0(file.path(outDir, species(thisEv)$id, id(thisEv)), '.csv'), row.names=FALSE)
        if(progress) {
            setTxtProgressBar(pb, value=e)
        }
    }
    if(length(naSp) > 0) {
        warning('No species id for events ', printN(naSp), call.=FALSE)
    }
    if(length(noDet) > 0) {
        warning('No detections found in events ', printN(naSp), call. = FALSE)
    }
    if(length(noBin) > 0) {
        warning('No binary data found for events ', printN(noBin), call.=FALSE)
    }
    if(progress) {
        cat('\n')
    }
    invisible(list(files=allFiles,
                   warnings=list(noSpecies=naSp,
                                 noDetections=noDet,
                                 noBinary=noBin)))
}

