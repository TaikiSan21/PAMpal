#' @title Set the Species Classification of an Acoustic Event
#'
#' @description Sets the \code{species} slot of an \linkS4class{AcousticEvent}
#'   object or list of objects
#'
#' @param x a \linkS4class{AcousticStudy} object, a list of \linkS4class{AcousticEvent}
#'   objects, or a single \linkS4class{AcousticEvent} object
#' @param type the type of classification to set, this is just a label within
#'   the \code{species} slot
#' @param method the method for assigning species to an event. Currently supports
#'   \code{pamguard}, which will use the 'eventType' or 'Text_Annotation' column
#'   to assign species, \code{manual} which will use \code{value} to assign
#'   species manually, or \code{reassign} which will use \code{value} to
#'   reassign an old species label to a new one
#' @param value optional argument required if \code{method} is set to manual or reassign.
#'   For \code{'manual'}, can either be a single value to assign to all events, or if assigning to
#'   a list a vector with length equal to the list. Can also be a dataframe
#'   with columns \code{event} and \code{species}, in which case species will
#'   be matched to corresponding event names instead of just relying on the
#'   order. If using this, please note the prefix OE or DGL present on most
#'   event numbers (see the \code{id} slot of your events).
#'   For \code{'reassign'}, \code{value} must be a data frame with columns
#'   \code{old} and \code{new}. Any events with species id in the \code{old} column
#'   of the dataframe will get reassigned to the corresponding id in the
#'   \code{new} column.
#' @return the same object as \code{x}, with species identifications assigned
#'   as an item named \code{type} in the \code{species} slot
#'
#' @author Taiki Sakai \email{taiki.sakai@@noaa.gov}
#'
#' @importFrom dplyr distinct
#' @importFrom stringr str_trim str_split
#' @importFrom RSQLite dbConnect dbDisconnect dbReadTable SQLite
#' @export
#'
setSpecies <- function(x, type='id', method=c('pamguard', 'manual', 'reassign'), value) {
    if(length(method) > 1) {
        cat('Please select a single assignment method.')
        methodPick <- menu(choices = method)
        if(methodPick == 0) {
            cat('No assignment method chosen, species cannot be assigned.')
            return(x)
        }
        method <- method[methodPick]
    }
    method <- match.arg(method[1], choices = c('pamguard', 'manual', 'am', 'reassign'))
    # wrote code for just events before study existed, just work on those
    if(is.AcousticStudy(x)) {
        acev <- events(x)
    # rest of code expects list for indexing because im bad and lazy
    } else if(is.AcousticEvent(x)) {
        acev <- list(x)
    # i know this is ugly dont judge me
    } else if(is.list(x)) {
        acev <- x
    }
    switch(method,
           'pamguard' = {
               spCol <- c('Text_Annotation', 'eventType', 'eventLabel')
               for(i in seq_along(acev)) {
                   sp <- sapply(detectors(acev[[i]]), function(y) {
                       hasCol <- spCol[spCol %in% colnames(y)]
                       unique(y[, hasCol])
                   })
                   sp <- unique(sp)
                   if(length(sp) > 1) {
                       spix <- menu(title = paste0('More than one species found for event ',
                                                   id(acev[[i]]), ', select one to assign:'),
                                    choices = sp)
                       if(spix == 0) {
                           warning('No species selected, assigning NA. Please fix later.')
                           sp <- NA_character_
                       } else {
                           sp <- sp[spix]
                       }
                   }
                   species(acev[[i]])[[type]] <- sp
               }
           },
           'manual' = {
               if(missing(value)) {
                   warning('Manual mode requires a "value" to set."')
                   return(x)
               }
               if(inherits(value, 'data.frame')) {
                   if(!all(c('species', 'event') %in% colnames(value))) {
                       warning('If "value" is a dataframe it must contain columns species and event.')
                       return(x)
                   }
                   allIds <- sapply(acev, id)
                   hasId <- allIds %in% value$event
                   if(!all(hasId)) {
                       warning('No match found for event(s) ',
                               paste0(allIds[!hasId], collapse=', '),
                               ' (Event names in "value" must match exactly)')
                   }
                   cat('Assigning species ids to ', sum(hasId), ' events.\n', sep='')
                   for(i in which(hasId)) {
                       species(acev[[i]])[[type]] <- value[value$event == id(acev[[i]]), 'species']
                   }
                # case when no a data frame
               } else {
                   if(length(value) != 1 &&
                      length(value) != length(acev)) {
                       warning('Length of "value" must be either 1 or the number of events.')
                       return(x)
                   }
                   if(length(value) == 1) {
                       value <- rep(value, length(acev))
                   }
                   for(i in seq_along(acev)) {
                       species(acev[[i]])[[type]] <- value[i]
                   }
               }
           },
           'am' = {
               specDf <- distinct(do.call(rbind, lapply(acev, function(oneAe) {
                   dbs <- files(oneAe)$db
                   events <- do.call(rbind, lapply(dbs, function(y) {
                       con <- dbConnect(y, drv=SQLite())
                       evs <- dbReadTable(con, 'Click_Detector_OfflineEvents')
                       dbDisconnect(con)
                       # browser()
                       evs <- evs[, c('UID', 'eventType', 'comment')]
                       evs$event <- paste0(gsub('\\.sqlite3', '', basename(y)),
                                              '.OE', as.character(evs$UID))
                       evs$eventType <- str_trim(evs$eventType)
                       evs$comment <- gsub('OFF EFF', '', evs$comment)
                       evs$comment <- gsub("[[:punct:]]", '', evs$comment)
                       evs$comment <- str_trim(evs$comment)
                       evs
                   }))
                   # events$event <- paste0('OE', as.character(events$UID))
                   events$species <- 'unid'
                   goodEvents <- c('BEAK', 'FORG')
                   events$species[events$eventType %in% goodEvents] <- str_split(events$comment[events$eventType %in% goodEvents],
                                                                         ' ', simplify=TRUE)[, 1]
                   events$species <- tolower(events$species)
                   events$species[events$species %in% c('mmme', 'mm')] <- 'unid'
                   events
               }
               )))
               specToAssign <- unique(specDf[specDf$event %in% sapply(acev, id), 'species'])
               if(length(specToAssign) > 0) {
                   cat('Assigning unique species: ', paste0(specToAssign, collapse = ', '), '.\n', sep = '')
               }
               acev <- setSpecies(acev, method = 'manual', type=type, value = specDf)
           },
           'reassign' = {
               if(missing(value)) {
                   warning('"reassign" mode requires a "value" dataframe.')
                   return(x)
               }
               colnames(value) <- tolower(colnames(value))
               if(!all(c('old', 'new') %in% colnames(value))) {
                   warning('Data frame must have columns "old" and "new" to reassign.')
                   return(x)
               }
               unchanged <- vector('character', length=0)
               for(i in seq_along(acev)) {
                   oldSpec <- species(acev[[i]])[[type]]
                   newSpec <- value[value$old == oldSpec, c('new')]
                   if(length(newSpec) == 0) {
                       unchanged <- ifelse(oldSpec %in% unchanged, unchanged, c(unchanged, oldSpec))
                       newSpec <- oldSpec
                   }
                   species(acev[[i]])[[type]] <- as.character(newSpec)
               }
               if(length(unchanged) > 0) {
                   cat(length(unchanged), ' species (', paste0(unchanged, collapse=', '), ') ',
                       'were not in reassignment dataframe, they have not been changed.', sep='')
               }
           },
           warning('Method ', method, ' not supported.')
    )
    if(is.AcousticStudy(x)) {
        events(x) <- acev
        return(x)
    }
    if(is.AcousticEvent(x)) {
        return(acev[[1]])
    }
    acev
}
