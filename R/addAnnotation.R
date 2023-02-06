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
#' @param x an \linkS4class{AcousticStudy} or \linkS4class{AcousticEvent} object
#' @param anno an annotation dataframe
#' @param verbose logical flag to print messages
#'
#' @return The same object as \code{x} with an \code{$annotation} item added
#'   the the \code{ancillary} slot of each event in \code{x}
#'
#' @author Taiki Sakai \email{taiki.sakai@@noaa.gov}
#'
#' @importFrom dplyr bind_rows arrange join_by
#'
#' @export
#'
addAnnotation <- function(x, anno, verbose=TRUE) {
    if(is.AcousticEvent(x)) {
        if(!id(x) %in% anno$event) {
            return(x)
        }
        thisAnno <- anno[anno$event == id(x), ]
        oldAnno <- ancillary(x)$annotation
        ancillary(x)$annotation <- safeListAdd(oldAnno, thisAnno)
        return(x)
    }
    if(!is.AcousticStudy(x)) {
        stop('"x" must be an AcousticStudy or AcousticEvent')
    }
    nAdded  <- 0
    for(e in seq_along(events(x))) {
        if(!id(x[[e]]) %in% anno$event) {
            next
        }
        x[[e]] <- addAnnotation(x[[e]], anno)
        nAdded <- nAdded + 1
    }
    if(verbose) {
        cat('Added annotations to ', nAdded, ' out of ',
            length(events(x)), ' events', sep='')
    }
    x
}

#' @param specMap data.frame to map species ids in \code{x} to names to be used
#'   for the annotation (ex. from \code{'ZC'} to \code{'Ziphius cavirostris'}.
#'   Dataframe must have columns \code{old} and \code{new}
#' @param mode one of \code{'event'} or \code{'detection'} to create annotation for
#'   events or detections
#' @param interactive logical flag to fill annotation data interactively
#'   (not recommended)
#' @param \dots additional named arguments to fill in annotation data. If names match
#'   a column in the annotation, that value will be used for all events or detections
#'   in the annotation
#'   
#' @export
#' @rdname addAnnotation
#'
prepAnnotation <- function(x, specMap=NULL, mode=c('event', 'detection'), interactive=FALSE, ...) {
    # Check first if GPS has been added
    dets <- getDetectorData(x)
    detNos <- sapply(dets, nrow)
    mode <- match.arg(mode)
    if(all(detNos == 0)) {
        warning('No detections found, no annotation added.')
        return(x)
    }
    dets <- dets[detNos > 0]
    if(!all(c('UTC', 'Longitude', 'Latitude') %in% colnames(dets[[1]]))) {
        # manually add lat/long?? can just set them to NA here
        stop('No GPS data found, please add first with "addGps".')
    }
    
    annoList <- getAnnoTemplate(notes=FALSE)
    
    fillEventAnno <- function(event, anno, interactive=interactive, mode) {
        coordsOnly <- bind_rows(lapply(detectors(event), function(d) {
            d[, c('UTC', 'Longitude', 'Latitude', 'UID')]
        }))
        
        switch(mode,
               'event' = {
                   useIx <- which.min(coordsOnly$UTC)[1]
                   anno$time_end <- max(coordsOnly$UTC)
               },
               'detection' = {
                   useIx <- 1:nrow(coordsOnly)
                   anno$time_end <- coordsOnly$UTC
               }
        )
        anno$lon <- coordsOnly$Longitude[useIx]
        anno$lat <- coordsOnly$Latitude[useIx]
        anno$time_start <- coordsOnly$UTC[useIx]
        
        thisSpec <- species(event)$id
        if(!is.null(specMap)) {
            thisSpec <- specMap$new[specMap$old == thisSpec]
        }
        anno$taxon <- thisSpec
        switch(mode,
               'detection' = {
                   freq <- getFrequency(event)
                   anno$freq_low <- freq$freqMin
                   anno$freq_high <- freq$freqMax
                   anno$matchId <- paste0(id(event), '.', coordsOnly$UID)
               },
               'event' ={
                   anno$freq_low <- NA
                   anno$freq_high <- NA
                   anno$matchId <- id(event)
               }
        )
        
        if(interactive) {
            cat(paste0('Adding data for event ', id(event), '...\n'))
            anno <- fillAnno(anno)
        }
        anno
    }
    
    allAnno <- lapply(events(x), function(e) {
        fillEventAnno(e, annoList, interactive, mode=mode)
    })
    allAnno <- bind_rows(allAnno, .id='event')
    
    dotArgs <- list(...)
    for(d in names(dotArgs)) {
        if(!d %in% names(allAnno)) {
            next
        }
        allAnno[[d]] <- dotArgs[[d]]
    }
    
    allAnno
}

#' @export
#' @rdname addAnnotation
#'
getAnnotation <- function(x) {
    if(is.AcousticEvent(x)) {
        return(ancillary(x)$annotation)
    }
    if(!is.AcousticStudy(x)) {
        stop('"x" must be an AcousticStudy or AcousticEvent.')
    }
    anno <- lapply(events(x), getAnnotation)
    bind_rows(anno, .id='event')
}

#' @export
#' @rdname addAnnotation
#' 
checkAnnotation <- function(x) {
    allCols <- names(getAnnoTemplate())
    missCols <- allCols[!allCols %in% names(x)]
    if(length(missCols) > 0) {
        warning(length(missCols), ' field(s) were not present: ',
                paste0(missCols, collapse=', '), call.=FALSE)
        for(c in missCols) {
            x[[c]] <- NA
        }
    }
    mandatory <- c('source', 'annotator', 'recording_url', 'contact')
    # non-mand
    notMand <- c('annotation_date', 'annotation_info_url', 'recording_info_url',
                 'time_start', 'time_end', 'freq_low', 'freq_high', 'species',
                 'type', 'lat', 'lon', 'alt', 'effort', 'additional')
    # notMand <- c('recording_info_url', 'annotation_info_url', 'type',
    #              'annotation_date', 'alt', 'freq_low','freq_high', 'effort',
    #              'additional')
    hasNa <- apply(x, 2, function(c) any(is.na(c)))
    mandMiss <- hasNa[mandatory]
    otherMiss <- hasNa[notMand]
    if(any(mandMiss)) {
        warning('Fill in data for ', sum(mandMiss), ' MANDATORY field(s): ',
                paste0(mandatory[mandMiss], collapse=', '), call.=FALSE)
    }
    if(any(otherMiss)) {
        message('Also missing data for ', sum(otherMiss), ' optional field(s): ',
                paste0(names(which(otherMiss)), collapse=', '), call.=FALSE)
    }
    negStart <- x$time_start < 0
    if(any(negStart)) {
        warning('Found ', sum(negStart, na.rm=TRUE), ' negative time_start values',
                call.=FALSE)
    }
    negEnd <- x$time_end < 0
    if(any(negEnd)) {
        warning('Found ', sum(negEnd, na.rm=TRUE), ' negative time_end values',
                call.=FALSE)
    }
    startAfter <- x$time_start > x$time_end
    if(any(startAfter)) {
        warning('Found ', sum(startAfter, na.rm=TRUE), ' occasions where time_start is ',
                'later than time_end', call.=FALSE)
    }
    freqLowHigher <- x$freq_low > x$freq_high
    if(!all(is.na(freqLowHigher)) &&
       sum(freqLowHigher, na.rm=TRUE) > 0) {
        warning('Found ', sum(freqLowHigher, na.rm=TRUE), ' occasions where freq_low is ',
                'higher than freq_high', call.=FALSE)
    }
    x
}

getAnnoTemplate <- function(notes=FALSE) {
    if(!notes) {
        return(
            list(
                source = NA, # github? --- study *
                source_id = NA, #unclear, unique --- event ?? or is this just a number identifying source *
                annotator = NA, # ask/fillin alg or person --- event repeat *
                annotation_id = NA, #unique for a source --- event repeat *
                annotation_date = NA, # ask/fillin
                annotation_info_url = NA, # github url?  --- event
                recording_url = NA, # github url --- event *
                recording_info_url = NA, # github folder? --- event
                time_start = TRUE, #getEventTime in utils reports named vector start/stop per event
                time_end = TRUE,
                freq_low = TRUE, # only if freq is in data we ave? could be optional i guess. "peak" in standardClick, "freqMin/Max" in rocca
                freq_high = TRUE,
                species = TRUE, # can ask for each species found in species id. was prev taxon?
                type = NA, # probably ask, could do whistle/click/bp based on detector
                lat = TRUE, # all trues are able to be autod, these do by start time
                lon = TRUE,
                alt = NA, # blank unless we have hydro depth, not currently looking for this when i read but could
                effort = NA, # lol
                contact = NA, # ask/fill in --- event repeat *
                additional = NA # ???
            )
        )
    } else {
        return(
            list(
                source = 'Where the annotated recording comes from (e.g. GitHub, bio.acousti.ca)', # github? --- study
                source_id = 'Identifier for the source recording, must be unique', #unclear, unique --- event ?? or is this just a number identifying source
                # ^ APPEARS TO MATCH TO UNIQUE RECORDING URL PER EDS EXAMPLES^
                annotator = 'Person or algorithm who made the annotations', # ask/fillin alg or person --- event repeat
                annotation_id = 'Identifier for the annotations, must be unique', #unique for a source --- event repeat
                annotation_date = 'Date/time the annotations were made', # ask/fillin
                annotation_info_url = 'URL of additional information relating to the annotation', # github url?  --- event
                recording_url = 'URL of the recording', # github url --- event
                recording_info_url = 'URL of additional information relating to the recording', # github folder? --- event
                # time_start = TRUE, #getEventTime in utils reports named vector start/stop per event.
                # time_end = TRUE, # These are SECONDS WITHIN RECORDING
                # freq_low = TRUE, # only if freq is in data we ave? could be optional i guess. "peak" in standardClick, "freqMin/Max" in rocca
                # freq_high = TRUE,
                # species = NA, # can ask for each species found in species id
                type = 'Type of call in the recording (e.g. click, whistle, noise)', # probably ask, could do whistle/click/bp based on detector
                # lat = TRUE, # all trues are able to be autod, these do by start time
                # lon = TRUE,
                alt = 'Recorder depth', # blank unless we have hydro depth, not currently looking for this when i read but could
                effort = 'Description of the effort used to annotate the file', # lol
                contact = 'Contact person for queries relating to the annotations (e-mail)', # ask/fill in --- event repeat
                additional = 'Any additional information' # ???
            )
        )
    }
}

# prepAnnotation <- function(x, interactive=TRUE) {
#     annoList <- getAnnoTemplate(notes=FALSE)
#     if(isFALSE(interactive)) {
#         return(annoList)
#     }
#     if(is.na(species(x[[1]])$id)) {
#         spStop <- menu(
#             title = paste0('Species have not been yet been assigned with "setSpecies"',
#                            ' would you like to continue?'),
#             choices = c('Yes, I will set "taxon" manually.',
#                         'No, stop preparing annotation.')
#         )
#         if(spStop %in% c(0, 2)) {
#             stop('Exiting prepAnnotation', call.=FALSE)
#         }
#         if(spStop == 1) {
#             annoList$taxon <- NA
#         }
#     }
#     notes <- getAnnoTemplate(notes=TRUE)
#     toFill <- sapply(annoList, is.na)
#     
#     for(i in names(annoList)[toFill]) {
#         thisChoice <- menu(
#             title = paste0('Will "', i, '" be the same for all events?',
#                            ' (', notes[[i]], ')'),
#             choices = c('Yes, same for all events.',
#                         'No, could be different for some events.',
#                         'This will not be filled in for ANY events.')
#         )
#         if(thisChoice %in% c(0, 2)) {
#             next
#         }
#         if(thisChoice == 1) {
#             cat(paste0('Enter value for "', i, '": ', notes[[i]]))
#             val <- readline()
#             if(val == '') {
#                 val <- 'NA'
#             }
#             annoList[[i]] <- val
#             next
#         }
#         if(thisChoice == 3) {
#             annoList[[i]] <- FALSE
#         }
#     }
#     if(isTRUE(annoList$taxon)) {
#         SPMATCH <- unique(sapply(events(x), function(e) species(e)$id))
#         names(SPMATCH) <- SPMATCH
#         cat(paste0('Found ', length(SPMATCH), ' species codes, provide the full names for annotation.\n'))
#         for(s in seq_along(SPMATCH)) {
#             cat(paste0('What should species ', SPMATCH[s], ' be called in the annotation?\n'))
#             SPMATCH[s] <- readline()
#         }
#         annoList$taxon <- SPMATCH
#     }
#     annoList
# }

getFrequency <- function(x) {
    if(!is.AcousticEvent(x)) {
        stop('Must be an AcousticEvent')
    }
    dets <- detectors(x)
    # freqNames <- c('peak', 'freqMin', 'freqMax')
    freqNames <- c('freqMin', 'freqMax')
    freqMin <- vector('numeric', length=0)
    freqMax <- vector('numeric', length=0)
    for(i in seq_along(dets)) {
        hasFreq <- freqNames %in% colnames(dets[[i]])
        if(!any(hasFreq)) {
            freqMin <- c(freqMin, rep(NA, nrow(dets[[i]])))
            freqMax <- c(freqMax, rep(NA, nrow(dets[[i]])))
            next
        }
        for(f in freqNames[hasFreq]) {
            thisFreq <- dets[[i]][[f]]
            # peak is in khz, others in hz
            # if(f == 'peak') {
            #     thisFreq <- thisFreq * 1e3
            # }
            if(f == 'freqMin') {
                freqMin <- c(freqMin, thisFreq)
            }
            if(f == 'freqMax') {
                freqMax <- c(freqMax, thisFreq)
            }
            # freq <- c(freq, thisFreq)
        }
    }
    if(length(freqMin) == 0) {
        freqMin <- NA
    }
    if(length(freqMax) == 0) {
        freqMax <- NA
    }
    list(freqMin=freqMin, freqMax=freqMax)
}

# if NA we ask to fill providing note
# if TRUE it fills auto
# used only for interactive fill
fillAnno <- function(anno=NULL, recheck=FALSE) {
    if(is.null(anno)) {
        anno <- getAnnoTemplate(notes=FALSE)
    }
    notes <- getAnnoTemplate(notes=TRUE)
    toFill <- sapply(anno, is.na)
    if(recheck) {
        # do something here, cant make all true cuz lat/long on autofill etc.
    }
    if(!any(toFill)) {
        return(anno)
    }
    for(i in names(anno)[toFill]) {
        cat(paste0('Enter value for "', i, '": ', notes[[i]]))
        val <- readline()
        if(val == '') {
            val <- 'NA'
        }
        anno[[i]] <- val
    }
    anno
}

# getAnnoId <- function(field=c('source_id', 'annotation_id'), match) {
#     return(0)
#     matchField <- switch(match.arg(field),
#                          'source_id' = '?recording_url=',
#                          'annotation_id' = '?source='
#     )
#     rj <- rjson::fromJSON(file=paste0('https://api.audioblast.org/annotations/',
#                                       matchField,
#                                       match))
#     ids <- sapply(rj, function(x) x[[field]])
#     ids <- suppressWarnings(as.numeric(ids))
#     ids <- ids[!is.na(ids)]
#     if(length(ids) == 0) {
#         ids <- return(0)
#     }
#     ids
# }
#' @param file file to write a CSV of the annotations to, if \code{NULL} (default)
#'   then no file will be written
#' @export
#' @rdname addAnnotation
#'
export_annomate <- function(x, file=NULL) {
    if(is.AcousticStudy(x)) {
        anno <- getAnnotation(x)
    } else {
        anno <- x
    }
    annomateCols <- names(getAnnoTemplate())
    anno <- anno[annomateCols]
    checkAnnotation(anno)
    if(!is.null(file)) {
        write.csv(x, file=file, row.names = FALSE)
    }
    anno
}

# originally developed for figshare uploads of data created with
# PAMpal::writeEventClips, so defaults match these patterns
# recording_url, filename, (matchId), (filestart)
#' @param rec dataframe of recording url information. Must have column
#'   \code{recording_url}. If clips were created using \link[PAMpal]{writeEventClips},
#'   then must have column \code{filename} containing the wav file names. Other
#'   column names will be automatically parsed from there. If wav files are from
#'   another source, must contain columns \code{matchId} 
#'
#' @export
#' @rdname addAnnotation
#' 
matchRecordingUrl <- function(anno, rec) {
    # matchId is either event name or event.UID for detections
    if(!'matchId' %in% names(rec)) {
        rec[['matchId']] <- parseEventClipName(rec[['filename']], 'event')
    }
    if(!'filestart' %in% names(rec)) {
        rec[['filestart']] <- parseEventClipName(rec[['filename']], 'time')
    }
    rec <- rec[c('matchId', 'recording_url', 'filestart')]
    origOrder <- names(anno)
    anno[['recording_url']] <- NULL
    # dont want to join on time if not necessary in case it doesnt work
    test <- left_join(anno, rec, by=join_by('matchId' == 'matchId'), multiple='all')
    if(nrow(test) > nrow(anno)) {
        anno <- left_join(anno, rec, by=join_by('matchId'=='matchId', closest('time_start' >= 'filestart')))
    } else {
        anno <- test
    }
    anno$time_start <- round(as.numeric(difftime(anno$time_start, anno$filestart, units='secs')), 3)
    anno$time_end <- round(as.numeric(difftime(anno$time_end, anno$filestart, units='secs')), 3)
    anno[origOrder]
}

globalVariables('closest')
