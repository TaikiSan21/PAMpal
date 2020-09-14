#' @title Add Annotation Data to an AcousticStudy Object
#'
#' @description Adds annotation data to an \linkS4class{AcousticStudy}
#'   object, usually in preparation for exporting to an external source.
#'   Most pieces of the annotation form will be filled in by the user when
#'   the function is called, but \code{lat}, \code{lon}, \code{time_start},
#'   \code{time_end}, \code{freq_low}, \code{freq_high}, \code{source_id},
#'   and \code{annotation_id} will be filled in automatically based on data
#'   in each \linkS4class{AcousticEvent}. Annotations are stored for each
#'   event in the \code{ancillary} slot.
#'
#' @param x an \linkS4class{AcousticStudy} object
#' @param annoList an annotation list created by \code{getAnnotation}.
#'   Default \code{NULL} will create a blank list
#'
#' @return The same object as \code{x} with an \code{$annotation} item added
#'   the the \code{ancillary} slot of each event in \code{x}
#'
#' @author Taiki Sakai \email{taiki.sakai@@noaa.gov}
#'
#' @importFrom dplyr bind_rows arrange
#' @importFrom rjson fromJSON
#' @importFrom rlang .data
#'
#' @export
#'
addAnnotation <- function(x, annoList=NULL) {
    # Check first if GPS has been added
    dets <- getDetectorData(x)
    detNos <- sapply(dets, nrow)
    if(all(detNos == 0)) {
        warning('No detections found, no annotation added.')
        return(x)
    }
    dets <- dets[detNos > 0]
    if(!all(c('UTC', 'Longitude', 'Latitude') %in% colnames(dets[[1]]))) {
        # manually add lat/long?? can just set them to NA here
        stop('No GPS data found, please add first with "addGps".')
    }

    SPMATCH <- unique(sapply(events(x), function(e) species(e)$id))
    names(SPMATCH) <- SPMATCH
    cat(paste0('Found ', length(SPMATCH), ' species codes, provide the full names for annotation.\n'))
    for(s in SPMATCH) {
        cat(paste0('What should species ', SPMATCH[[s]], ' be called in the annotation?\n'))
        SPMATCH[[s]] <- readline()
    }
    if(is.null(annoList)) {
        annoList <- getAnnotation(x)
    }
    # Fill common stuff first
    common <- c('source', 'annotator', 'contact')
    cat('Adding data common to all annotations...\n')
    annoList[common] <- fillAnno(annoList[common])
    # then pull IDs from annomate API. These will be 0 if none present, so next is 1.
    SOURCEID <- getAnnoId(field='source_id', match=annoList['recording_url'])
    ANNOTATIONID <- getAnnoId(field='annotation_id', match=annoList['source'])
    # this is event level stuff, do later!

    fillEventAnno <- function(event, annoList) {
        coordsOnly <- bind_rows(lapply(detectors(event), function(d) {
            d[, c('UTC', 'Longitude', 'Latitude')]
        }))
        coordsOnly <- arrange(coordsOnly, .data$UTC)
        annoList$lon <- coordsOnly$Longitude[1]
        annoList$lat <- coordsOnly$Latitude[1]
        annoList$time_start <- coordsOnly$UTC[1]
        annoList$time_end <- coordsOnly$UTC[nrow(coordsOnly)]
        nextSourceId <- max(SOURCEID) + 1
        annoList$source_id <- nextSourceId
        SOURCEID <<- c(SOURCEID, nextSourceId)
        nextAnnoId <- max(ANNOTATIONID) + 1
        annoList$annotation_id <- nextAnnoId
        ANNOTATIONID <<- c(ANNOTATIONID, nextAnnoId)
        annoList$taxon <- SPMATCH[[species(event)$id]]
        freq <- getFrequency(event)
        annoList$freq_low <- min(freq)
        annoList$freq_high <- max(freq)
        cat(paste0('Adding data for event ', id(event), '...\n'))
        annoList <- fillAnno(annoList)
        event@ancillary$annotation <- annoList
        event
    }

    for(e in seq_along(events(x))) {
        events(x)[[e]] <- fillEventAnno(events(x)[[e]], annoList)
    }
    x
    # on export, have it double check each field it has already found, have override option for programattc
    # have option to add or not annotation when processing data
}

#' @param notes logical flag whether to pull a version of the list that only contains notes,
#'   should only be \code{TRUE} for internal use
#' @export
#' @rdname addAnnotation
#'
getAnnotation <- function(x=NULL, notes=FALSE) {
    if(is.AcousticEvent(x)) {
        anno <- x@ancillary$annotation
        if(is.null(anno)) {
            return(getAnnotation(x=NULL, notes=FALSE))
        }
        return(anno)
    }
    if(!notes) {
        return(
            list(
                source = NA, # github? --- study
                source_id = NA, #unclear, unique --- event ?? or is this just a number identifying source
                recording_url = NA, # github url --- event
                recording_info_url = NA, # github folder? --- event
                annotator = NA, # ask/fillin alg or person --- event repeat
                annotation_id = NA, #unique for a source --- event repeat
                contact = NA, # ask/fill in --- event repeat
                annotation_info_url = NA, # github url?  --- event
                taxon = NA, # can ask for each species found in species id
                type = NA, # probably ask, could do whistle/click/bp based on detector
                annotation_date = NA, # ask/fillin
                lat = TRUE, # all trues are able to be autod, these do by start time
                lon = TRUE,
                alt = NA, # blank unless we have hydro depth, not currently looking for this when i read but could
                time_start = TRUE, #getEventTime in utils reports named vector start/stop per event
                time_end = TRUE,
                freq_low = TRUE, # only if freq is in data we ave? could be optional i guess. "peak" in standardClick, "freqMin/Max" in rocca
                freq_high = TRUE,
                effort = NA, # lol
                additional = NA # ???
            )
        )
    } else {
        return(
            list(
                source = 'Where the annotated recording comes from (e.g. GitHub, bio.acousti.ca)', # github? --- study
                source_id = 'Identifier for the source recording, must be unique', #unclear, unique --- event ?? or is this just a number identifying source
                # ^ APPEARS TO MATCH TO UNIQUE RECORDING URL PER EDS EXAMPLES^
                recording_url = 'URL of the recording', # github url --- event
                recording_info_url = 'URL of additional information relating to the recording', # github folder? --- event
                annotator = 'Person or algorithm who made the annotations', # ask/fillin alg or person --- event repeat
                annotation_id = 'Identifier for the annotations, must be unique', #unique for a source --- event repeat
                contact = 'Contact person for queries relating to the annotations (e-mail)', # ask/fill in --- event repeat
                annotation_info_url = 'URL of additional information relating to the annotation', # github url?  --- event
                # taxon = NA, # can ask for each species found in species id
                type = 'Type of call in the recording (e.g. click, whistle, noise)', # probably ask, could do whistle/click/bp based on detector
                annotation_date = 'Date/time the annotations were made', # ask/fillin
                # lat = TRUE, # all trues are able to be autod, these do by start time
                # lon = TRUE,
                alt = 'Recorder depth', # blank unless we have hydro depth, not currently looking for this when i read but could
                # time_start = TRUE, #getEventTime in utils reports named vector start/stop per event.
                # time_end = TRUE, # These are SECONDS WITHIN RECORDING
                # freq_low = TRUE, # only if freq is in data we ave? could be optional i guess. "peak" in standardClick, "freqMin/Max" in rocca
                # freq_high = TRUE,
                effort = 'Description of the effort used to annotate the file', # lol
                additional = 'Any additional information' # ???
            )
        )
    }
}

getFrequency <- function(x) {
    if(!is.AcousticEvent(x)) {
        stop('Must be an AcousticEvent')
    }
    dets <- detectors(x)
    freqNames <- c('peak', 'freqMin', 'freqMax')
    freq <- vector('numeric', length=0)
    for(i in seq_along(dets)) {
        hasFreq <- freqNames %in% colnames(dets[[i]])
        if(!any(hasFreq)) next
        for(f in freqNames[hasFreq]) {
            thisFreq <- dets[[i]][[f]]
            # peak is in khz, others in hz
            if(f == 'peak') {
                thisFreq <- thisFreq * 1e3
            }
            freq <- c(freq, thisFreq)
        }
    }
    if(length(freq) == 0) {
        return(NA)
    }
    freq
}

fillAnno <- function(annoList=NULL, recheck=FALSE) {
    if(is.null(annoList)) {
        annoList <- getAnnotation(notes=FALSE)
    }
    notes <- getAnnotation(notes=TRUE)
    toFill <- sapply(annoList, is.na)
    if(recheck) {
        # do something here, cant make all true cuz lat/long on autofill etc.
    }
    if(!any(toFill)) {
        return(annoList)
    }
    for(i in names(annoList)[toFill]) {
        cat(paste0('Enter value for "', i, '": ', notes[[i]]))
        val <- readline()
        if(val == '') {
            val <- 'NA'
        }
        annoList[[i]] <- val
    }
    annoList
}

getAnnoId <- function(field=c('source_id', 'annotation_id'), match) {
    matchField <- switch(match.arg(field),
                         'source_id' = '?recording_url=',
                         'annotation_id' = '?source='
    )
    rj <- fromJSON(file=paste0('https://api.audioblast.org/annotations/',
                                      matchField,
                                      match))
    ids <- sapply(rj, function(x) x[[field]])
    ids <- suppressWarnings(as.numeric(ids))
    ids <- ids[!is.na(ids)]
    if(length(ids) == 0) {
       ids <- return(0)
    }
    ids
}
