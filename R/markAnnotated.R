#' @title Mark Detections as Annotated
#'
#' @description Marks detections within an \linkS4class{AcousticStudy} as
#'   being within the bounds of an annotation box. Annotations can either
#'   be read in from the "Spectrogram Annotation" module of PAMguard, or
#'   supplied as a separate dataframe. Detections must be entirely contained
#'   within the annotation bounds.
#'
#' @param x an AcousticStudy object
#' @param anno annotations to read from. If \code{NULL}, will be read in from
#'   the PAMguard database. If a data.frame, must have columns \code{start} and
#'   \code{end} in UTC, and column \code{id}. Can additionally have columns \code{fmin} and \code{fmax}
#'   to apply frequency bounds (values in Hz).
#' @param tBuffer additional buffer value to add on to annotation time bounds in
#'   seconds. If a single number, the number of seconds to extend the bounds by
#'   on the start and end of each annotation. Can also be a vector of two to
#'   extend different values on the start and end. This can be useful if original
#'   bounding boxes were drawn very close to the desired detections since any small
#'   portion of a signal outside the box will cause it to be excluded.
#' @param fBuffer additional buffer value to add to annotation frequency bounds in
#'   Hz. If a single number, the number of Hz to extend bounds by on lower and upper
#'   end of boxes. Can also be a vector of two to extend different values on lower and
#'   upper bounds. This can be useful if original bounding boxes were drawn very close
#'   to the desired detections since any small portion of a signal outside the box will
#'   cause it to be excluded.
#' @param table if \code{anno} is \code{NULL}, the name of the "Spectrogram Annotation"
#'   module table within the database.
#'
#' @details This adds new columns \code{inAnno} and \code{annoId} to all detector
#'   dataframes within the AcousticStudy. \code{inAnno} is a logical flag whether or not
#'   a given detection was fully contained in any annotation bounding box, and \code{annoId}
#'   lists the IDs of the boxes it matched. A detection is considered within an annotation
#'   only if it is entirely within the time and frequency bounds of the annotation. For
#'   GPL and whistle detections, the min and max frequency values are used. For click detections,
#'   only the peak frequency is used. For cesptrum detections, frequency bounds are ignored.
#'
#' @return the same object as \code{x}, but detectors have additional columns added
#'
#' @examples
#' data(exStudy)
#' annotation <- data.frame(start = min(getWhistleData(exStudy)$UTC),
#'                          fmin = c(16000, 17000),
#'                          fmax = c(17000, 18000))
#' annotation$end <- annotation$star + 1
#' exStudy <- markAnnotated(exStudy, annotation)
#' getWhistleData(exStudy)[c('UTC', 'duration', 'freqMin', 'freqMax', 'inAnno', 'annoId')]
#'
#' @author Taiki Sakai \email{taiki.sakai@@noaa.gov}
#'
#' @importFrom PAMmisc readSpecAnno
#'
#' @export
#'
markAnnotated <- function(x, anno=NULL, tBuffer=0, fBuffer=0, table='Spectrogram_Annotation') {
    # format buffers
    if(length(tBuffer) == 1) {
        tBuffer <- rep(tBuffer, 2) * c(-1, 1)
    }
    if(length(fBuffer) == 1) {
        fBuffer <- rep(fBuffer, 2) * c(-1, 1)
    }
    tBuffer <- abs(tBuffer) * c(-1, 1)
    fBuffer <- abs(fBuffer) * c(-1, 1)
    # format annotation list for different DB possiblities
    # and apply buffers
    allDb <- files(x)$db
    noAnnoTable <- character(0)
    if(is.null(anno)) {
        annoList <- lapply(allDb, function(d) {
            thisAnn <- readSpecAnno(d, table=table)
            if(is.null(thisAnn)) {
                noAnnoTable <<- c(noAnnoTable, basename(d))
                return(thisAnn)
            }
            thisAnn$start <- thisAnn$start + tBuffer[1]
            thisAnn$end <- thisAnn$end + tBuffer[2]
            thisAnn$start <- as.numeric(thisAnn$start)
            thisAnn$end <- as.numeric(thisAnn$end)
            thisAnn$fmin <- thisAnn$fmin + fBuffer[1]
            thisAnn$fmax <- thisAnn$fmax + fBuffer[2]
            thisAnn
        })
        names(annoList) <- basename(allDb)
    } else if(is.data.frame(anno)) {
        needCols <- c('start', 'end')
        if(!all(needCols %in% colnames(anno))) {
            stop('Annotations must have columns "start" and "end"')
        }
        if(!inherits(anno$start, 'POSIXct')) {
            anno$start <- parseUTC(anno$start)
        }
        if(!inherits(anno$end, 'POSIXct')) {
            anno$end <- parseUTC(anno$end)
        }
        if(!'fmin' %in% colnames(anno)) {
            anno$fmin <- 0
        }
        if(!'fmax' %in% colnames(anno)) {
            anno$fmax <- Inf
        }
        if(!'id' %in% colnames(anno)) {
            anno$id <- as.character(1:nrow(anno))
        }
        anno$start <- anno$start + tBuffer[1]
        anno$end <- anno$end + tBuffer[2]
        anno$start <- as.numeric(anno$start)
        anno$end <- as.numeric(anno$end)
        anno$fmin <- anno$fmin + fBuffer[1]
        anno$fmax <- anno$fmax + fBuffer[2]
        if('db' %in% colnames(anno)) {
            annoList <- split(anno, basename(anno$db))
            noAnnoTable <- c(noAnnoTable, setdiff(basename(allDb), names(annoList)))
        } else {
            annoList <- vector('list', length=length(allDb))
            for(a in seq_along(annoList)) {
                annoList[[a]] <- anno
            }
            names(annoList) <- basename(allDb)
        }
    }
    events(x) <- lapply(events(x), function(e) {
        thisAnn <- annoList[[basename(files(e)$db)]]
        for(d in seq_along(detectors(e))) {
            e[[d]] <- markOneDetector(e[[d]], thisAnn, type=attr(e[[d]], 'calltype'), event=id(e))
        }
        e
    })
    x <- .addPamWarning(x)
    x
}

markOneDetector <- function(x, sa, type=c('click', 'whistle', 'gpl', 'cepstrum'), event=NULL) {
    type <- match.arg(type)
    x$annoId <- ''
    x$inAnno <- FALSE
    if(is.null(sa)) {
        return(x)
    }
    x$numTime <- as.numeric(x$UTC)
    naUID <- character(0)
    for(i in 1:nrow(sa)) {
        inThis <- switch(
            type,
            'click' = {
                x$numTime >= sa$start[i] &
                    x$numTime <= sa$end[i] &
                # x$UTC >= sa$start[i] &
                    # x$UTC <= sa$end[i] &
                    x$peak * 1e3 <= sa$fmax[i] &
                    x$peak * 1e3 >= sa$fmin[i]
            },
            'whistle' = {
                # x$UTC >= sa$start[i] &
                    # (x$UTC + x$duration) <= sa$end[i] &
                x$numTime >= sa$start[i] &
                    (x$numTime + x$duration) <= sa$end[i] &
                    x$freqMax <= sa$fmax[i] &
                    x$freqMin >= sa$fmin[i]
            },
            'gpl' = {
                # x$UTC >= sa$start[i] &
                #     (x$UTC + x$duration) <= sa$end[i] &
                x$numTime >= sa$start[i] &
                    (x$numTime + x$duration) <= sa$end[i] &
                    x$freqMax <= sa$fmax[i] &
                    x$freqMin >= sa$fmin[i]
            },
            'cepstrum' = {
                # x$UTC >= sa$start[i] &
                #     (x$UTC + x$duration) <= sa$end[i]
                x$numTime >= sa$start[i] &
                    (x$numTime + x$duration) <= sa$end[i]
            }
        )
        naComp <- is.na(inThis)
        naUID <- unique(c(naUID, x$UID[naComp]))
        inThis[naComp] <- FALSE
        # store ID of spec anno box for matches just in case we want it later
        x$annoId[inThis] <- paste0(x$annoId[inThis], sa$id[i], ',')
        # change any detections in this box to TRUE
        x$inAnno <- x$inAnno | inThis
    }
    if(length(naUID) > 0) {
        pamWarning(length(naUID), ' detections in event ', event, ' had NA values and could not be properly',
                   ' compared to annotations (UIDs ', naUID, ')')
    }
    x$numTime <- NULL
    # switch(type,
    #        'click' = {
    #            for(i in 1:nrow(sa)) {
    #                # check every detection to see if it is fully contained within bounding box 'i'
    #                inThis <- x$UTC >= sa$start[i] &
    #                    x$UTC <= sa$end[i] &
    #                    x$peak * 1e3 <= sa$fmax[i] &
    #                    x$peak * 1e3 >= sa$fmin[i]
    #                # store ID of spec anno box for matches just in case we want it later
    #                x$annoId[inThis] <- paste0(x$annoId[inThis], sa$id[i], ',')
    #                # change any detections in this box to TRUE
    #                x$inAnno <- x$inAnno | inThis
    #            }
    #        },
    #        'whistle' = {
    #            for(i in 1:nrow(sa)) {
    #                # check every detection to see if it is fully contained within bounding box 'i'
    #                inThis <- x$UTC >= sa$start[i] &
    #                    (x$UTC + x$duration) <= sa$end[i] &
    #                    x$freqMax <= sa$fmax[i] &
    #                    x$freqMin >= sa$fmin[i]
    #                # store ID of spec anno box for matches just in case we want it later
    #                x$annoId[inThis] <- paste0(x$annoId[inThis], sa$id[i], ',')
    #                # change any detections in this box to TRUE
    #                x$inAnno <- x$inAnno | inThis
    #            }
    #        },
    #        'gpl' = {
    #            for(i in 1:nrow(sa)) {
    #                # check every detection to see if it is fully contained within bounding box 'i'
    #                inThis <- x$UTC >= sa$start[i] &
    #                    (x$UTC + x$duration) <= sa$end[i] &
    #                    x$freqMax <= sa$fmax[i] &
    #                    x$freqMin >= sa$fmin[i]
    #                # store ID of spec anno box for matches just in case we want it later
    #                # if(anyNA(inThis)) browser()
    #                x$annoId[inThis] <- paste0(x$annoId[inThis], sa$id[i], ',')
    #                # change any detections in this box to TRUE
    #                x$inAnno <- x$inAnno | inThis
    #            }
    #        },
    #        'cepstrum' = {
    #            for(i in 1:nrow(sa)) {
    #                # check every detection to see if it is fully contained within bounding box 'i'
    #                inThis <- x$UTC >= sa$start[i] &
    #                    (x$UTC + x$duration) <= sa$end[i]
    #                # store ID of spec anno box for matches just in case we want it later
    #                x$annoId[inThis] <- paste0(x$annoId[inThis], sa$id[i], ',')
    #                # change any detections in this box to TRUE
    #                x$inAnno <- x$inAnno | inThis
    #            }
    #        }
    # )
    x$annoId <- gsub(',$', '', x$annoId)
    x
}
