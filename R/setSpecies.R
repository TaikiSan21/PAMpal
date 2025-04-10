#' @title Set the Species Classification of Events
#'
#' @description Sets the \code{species} slot of \linkS4class{AcousticEvent}
#'   objects within an \linkS4class{AcousticStudy}
#'
#' @param x a \linkS4class{AcousticStudy} object, a list of \linkS4class{AcousticEvent}
#'   objects, or a single \linkS4class{AcousticEvent} object
#' @param method the method for assigning species to an event. Currently supports
#'   \code{pamguard}, which will use the 'eventType' or 'Text_Annotation' column
#'   to assign species, \code{manual} which will use \code{value} to assign
#'   species manually, or \code{reassign} which will use \code{value} to
#'   reassign an old species label to a new one
#' @param value required only if \code{method} is set to 'manual' or 'reassign'.
#'   For \code{'manual'}, can either be a single value to assign to all events, or a
#'   vector with length equal to the number of events. Can also be a dataframe
#'   with columns \code{event} and \code{species}, in which case species will
#'   be matched to corresponding event names instead of just relying on the
#'   order. If using this, please note the prefix OE or DGL present on most
#'   event numbers (see the \code{id} slot of your events, or \code{names(events(x))}).
#'   For \code{'reassign'}, \code{value} must be a data frame with columns
#'   \code{old} and \code{new}. Any events with species id in the \code{old} column
#'   of the dataframe will get reassigned to the corresponding id in the
#'   \code{new} column.
#' @param type the type of classification to set, this is just a label within
#'   the \code{species} slot. Default \code{'id'} should typically not be changed
#'   since this is used by other functions
#'
#' @return the same object as \code{x}, with species identifications assigned
#'   as an item named \code{type} in the \code{species} slot
#'
#' @author Taiki Sakai \email{taiki.sakai@@noaa.gov}
#'
#' @examples
#'
#' # setting up example data
#' exPps <- new('PAMpalSettings')
#' exPps <- addDatabase(exPps, system.file('extdata', 'Example.sqlite3', package='PAMpal'))
#' exPps <- addBinaries(exPps, system.file('extdata', 'Binaries', package='PAMpal'))
#' exClick <- function(data) {
#'     standardClickCalcs(data, calibration=NULL, filterfrom_khz = 0)
#' }
#' exPps <- addFunction(exPps, exClick, module = 'ClickDetector')
#' exPps <- addFunction(exPps, roccaWhistleCalcs, module='WhistlesMoans')
#' exPps <- addFunction(exPps, standardCepstrumCalcs, module = 'Cepstrum')
#' exData <- processPgDetections(exPps, mode='db')
#' exData <- setSpecies(exData, method='pamguard')
#' species(exData)
#' exData <- setSpecies(exData, method='manual', value = c('sp1', 'sp2'))
#' species(exData)
#' exData <- setSpecies(exData, method='reassign',
#'                      value = data.frame(old='sp1', new='sp3'))
#' species(exData)
#'
#' @importFrom dplyr distinct
#' @importFrom RSQLite dbConnect dbDisconnect dbReadTable SQLite
#' @export
#'
setSpecies <- function(x, method=c('pamguard', 'manual', 'reassign'), value, type='id') {
    if(length(method) > 1) {
        methodPick <- menu(title='Please select a single assignment method.', choices = method)
        if(methodPick == 0) {
            pamWarning('No assignment method chosen, species cannot be assigned.')
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
                           pamWarning('No species selected, assigning NA. Please fix later.')
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
                   pamWarning('Manual mode requires a "value" to set."')
                   return(x)
               }
               if(inherits(value, 'data.frame')) {
                   # jank to make it work with eventId bc I have regrets
                   if('eventId' %in% colnames(value) &&
                      !'event' %in% colnames(value)) {
                       colnames(value)[colnames(value) == 'eventId'] <- 'event'
                   }
                   if(!all(c('species', 'event') %in% colnames(value))) {
                       pamWarning('If "value" is a dataframe it must contain columns "species" and "event" or "eventId".')
                       return(x)
                   }
                   allIds <- sapply(acev, id)
                   hasId <- allIds %in% value$event
                   if(!all(hasId)) {
                       message('No match found for event(s) ',
                               printN(allIds[!hasId], 6),
                               ' (Event names in "value" must match exactly)')
                   }
                   # cat('Assigning species ids to ', sum(hasId), ' events.\n', sep='')
                   for(i in which(hasId)) {
                       # species(acev[[i]])[[type]] <- value[value$event == id(acev[[i]]), 'species', drop=TRUE]
                       species(acev[[i]])[[type]] <- value[['species']][value$event == id(acev[[i]])]
                   }
                # case when no a data frame
               } else {
                   if(length(value) != 1 &&
                      length(value) != length(acev)) {
                       pamWarning('Length of "value" must be either 1 or the number of events.')
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
               specDf <- bind_rows(lapply(acev, function(oneAe) {
                   list(event = id(oneAe),
                                  eventType = strsplitboth(getClickData(oneAe)$eventLabel[1]),
                                  comment = ancillary(oneAe)$eventComment)
                   # dbs <- files(oneAe)$db
                   # events <- do.call(rbind, lapply(dbs, function(y) {
                   #     con <- dbConnect(y, drv=SQLite())
                   #     evs <- dbReadTable(con, 'Click_Detector_OfflineEvents')
                   #     dbDisconnect(con)
                   #     # browser()
                   #     evs <- evs[, c('Id', 'eventType', 'comment')]
                   #     evs$event <- paste0(gsub('\\.sqlite3', '', basename(y)),
                   #                            '.OE', as.character(evs$Id))
                   #     evs$eventType <- strsplitboth(evs$eventType)
                   #     evs$comment <- gsub('OFF EFF', '', evs$comment)
                   #     evs$comment <- gsub("[[:punct:]]", '', evs$comment)
                   #     evs$comment <- strsplitboth(evs$comment)
                   #     evs
                   # }))
                   # events$event <- paste0('OE', as.character(events$UID))
                   # events$species <- 'unid'
                   # events$species[events$eventType %in% goodEvents] <- str_split(events$comment[events$eventType %in% goodEvents],
                   #                                                       ' ', simplify=TRUE)[, 1]
                   # events$species <- tolower(events$species)
                   # events$species[events$species %in% c('mmme', 'mm')] <- 'unid'
                   # events
               }))
               specDf$comment <- gsub('OFF EFF', '', specDf$comment)
               specDf$comment <- gsub("[[:punct:]]", '', specDf$comment)
               specDf$comment <- strsplitboth(specDf$comment)
               specDf$species <- 'unid'
               goodEvents <- c('BEAK', 'FORG')
               specDf$species <- 'unid'
               # specDf$species[specDf$eventType %in% goodEvents] <- str_split(specDf$comment[specDf$eventType %in% goodEvents],
               #                                                               ' ', simplify=TRUE)[, 1]
               specDf$species[specDf$eventType %in% goodEvents] <- sapply(specDf$comment[specDf$eventType %in% goodEvents], 
                                                                          function(x) {
                                                                              strsplit(x, ' ')[[1]][1]
                                                                          }, USE.NAMES = FALSE)
                   
               specDf$species <- tolower(specDf$species)
               specDf$species[specDf$species %in% c('mmme', 'mm')] <- 'unid'

               # specToAssign <- unique(specDf[specDf$event %in% sapply(acev, id), 'species'])
               specToAssign <- unique(specDf$species)
               if(length(specToAssign) > 0) {
                   cat('Assigning unique species: ', paste0(specToAssign, collapse = ', '), '.\n', sep = '')
               }
               acev <- setSpecies(acev, method = 'manual', type=type, value = specDf)
           },
           'reassign' = {
               if(missing(value)) {
                   pamWarning('"reassign" mode requires a "value" dataframe.')
                   return(x)
               }
               colnames(value) <- tolower(colnames(value))
               if(!all(c('old', 'new') %in% colnames(value))) {
                   pamWarning('Data frame must have columns "old" and "new" to reassign.')
                   return(x)
               }
               unchanged <- vector('character', length=0)
               for(i in seq_along(acev)) {
                   oldSpec <- species(acev[[i]])[[type]]
                   if(is.na(oldSpec)) {
                       if(any(is.na(value$old))) {
                           newSpec <- value$new[is.na(value$old)[1]]
                       } else {
                           newSpec <- NA
                       }
                   } else {
                       newSpec <- value[value$old == oldSpec, c('new')]
                   }
                   if(length(newSpec) == 0) {
                       unchanged <- c(unchanged, oldSpec)
                       newSpec <- oldSpec
                   }
                   species(acev[[i]])[[type]] <- as.character(newSpec)
               }
               unchanged <- unique(unchanged)
               if(length(unchanged) > 0) {
                   message(length(unchanged), ' species (', printN(unchanged, 6), ') ',
                       'were not in reassignment dataframe, they have not been changed.', sep='')
               }
           },
           pamWarning('Method ', method, ' not supported.')
    )
    if(is.AcousticStudy(x)) {
        events(x) <- acev
        x <- .addPamWarning(x)
        return(x)
    }
    if(is.AcousticEvent(x)) {
        return(acev[[1]])
    }
    acev
}
