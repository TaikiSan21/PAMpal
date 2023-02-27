#obis funs
library(sf)
library(dplyr)
library(PAMpal)
library(lubridate)
library(worrms)

# converts GPS dataframe to spatial linestring, only taking a point
# every "interval" days
gpsToLinestring <- function(x, interval=1) {
    x <- arrange(x, UTC)
    x$dayDiff <- as.numeric(difftime(x$UTC, x$UTC[1], units='days'))
    x$dayDiff <- floor(x$dayDiff / interval)
    less <- bind_rows(lapply(split(x, x$dayDiff), function(d) {
        d[1, ]
    }))
    st_linestring(matrix(c(less$Longitude, less$Latitude), ncol=2))
}

# return just some summary data for each event in acousticstudy
eventSummary <- function(x, measures=TRUE) {
    evData <- lapply(events(x), function(e) {
        result <- list()
        result$eventId <- id(e)
        coordsOnly <- bind_rows(lapply(detectors(e), function(d) {
            select(d, any_of(c('UTC', 'Longitude', 'Latitude')))
        }))
        coordsOnly <- arrange(coordsOnly, 'UTC')
        result$start <- coordsOnly$UTC[1]
        result$end <- coordsOnly$UTC[nrow(coordsOnly)]
        if(all(c('Latitude', 'Longitude') %in% names(coordsOnly))) {
            result$Latitude <- coordsOnly$Latitude[1]
            result$Longitude <- coordsOnly$Longitude[1]
        }
        result$species <- species(e)$id
        result$nClicks <- nClicks(e, distinct = TRUE)
        result$nWhistles <- nWhistles(e)
        result$nGPL <- nGPL(e)
        result$nCepstrum <- nCepstrum(e)
        result$db <- files(e)$db
        if(measures) {
            meas <- getMeasures(e)
            for(m in names(meas)) {
                if(m == 'id') {
                    next
                }
                result[[m]] <- meas[[m]]
            }
        }
        result
    })
    bind_rows(evData)
}

# create obis event dataframe from PAMpal
formatObisEvent <- function(gps, dbPattern, ...) {
    result <- bind_rows(lapply(split(gps, gps$db), function(d) {
        start <- min(d$UTC)
        end <- max(d$UTC)
        ls <- gpsToLinestring(d, interval=1)
        centroid <- st_centroid(ls)
        wkt <- st_as_text(ls)
        deployment <- gsub(dbPattern[1], dbPattern[2], basename(d$db[1]))
        start8601 <- paste0(lubridate::format_ISO8601(start), 'Z')
        end8601 <- paste0(lubridate::format_ISO8601(end), 'Z')
        id <- paste0(unique(deployment), '_', start8601)
        # id <- deployment
        list(eventID = id,
             eventDate = paste0(start8601, '/', end8601),
             decimalLatitude = centroid[2],
             decimalLongitude = centroid[1],
             countryCode = 'US',
             footprintWKT = wkt,
             footprintSRS = 'epsg:4326',
             geodeticDatum = 'WGS84',
             georeferenceProtocol = NA,
             samplingProtocol = NA,
             db = basename(d$db[1]))
    }))
    dotArgs <- list(...)
    if(is.null(names(dotArgs))) {
        dotArgs <- unlist(dotArgs, recursive = FALSE)
    }
    for(d in names(dotArgs)) {
        if(!d %in% names(result)) {
            next
        }
        result[[d]] <- dotArgs[[d]]
    }
    result
}

# create obis occurrence re ords from PAMpal
formatObisOccur <- function(eventSummary, speciesMap, dbPattern, ...) {
    wormsRecords <- bind_rows(wm_record(unique(speciesMap$aphia)))
    wr <- select(wormsRecords, 
                 scientificName = scientificname, 
                 scientificNameID = lsid, 
                 kingdom, 
                 taxonRank = rank,
                 AphiaID)
    speciesMap <- left_join(speciesMap, wr, by=c('aphia'='AphiaID'))
    eventSummary <- eventSummary[eventSummary$species %in% speciesMap$ours, ]
    eventSummary$obisId <- gsub(dbPattern[1], dbPattern[2], basename(eventSummary$db))
    obOc <- data.frame(
        occurrenceID = paste0(eventSummary$obisId, '_', 
                              lubridate::format_ISO8601(eventSummary$start), '_',
                              gsub(' ', '_', eventSummary$species)),
        basisOfRecord = 'MachineObservation',
        eventDate = paste0(lubridate::format_ISO8601(eventSummary$start), 'Z/',
                           lubridate::format_ISO8601(eventSummary$end), 'Z'),
        decimalLatitude = eventSummary$Latitude,
        decimalLongitude = eventSummary$Longitude,
        occurrenceStatus = 'present',
        countryCode = 'US',
        geodeticDatum = 'WGS84',
        # coordinateUncertaintyInMeters = 3130,
        institutionCode = NA,
        georeferenceProtocol = NA,
        samplingProtocol = NA,
        species = eventSummary$species,
        db = basename(eventSummary$db)
    )
    obOc <- left_join(obOc, speciesMap[c('ours', 'scientificName', 'scientificNameID',
                                         'kingdom', 'taxonRank', 'coordinateUncertaintyInMeters')],
                      by=join_by(species == ours))
    obOc$species <- NULL
    dotArgs <- list(...)
    if(is.null(names(dotArgs))) {
        dotArgs <- unlist(dotArgs, recursive = FALSE)
    }
    for(d in names(dotArgs)) {
        if(!d %in% names(obOc)) {
            next
        }
        obOc[[d]] <- dotArgs[[d]]
    }
    obOc
}

# do all event and occurrence steps directly from acoustic study
# and match their IDs appropriately
export_obis <- function(x, speciesMap, dbPattern=NULL, ...) {
    if(is.null(dbPattern) ||
       length(dbPattern) == 0) {
        dbPattern <- c('', '')
    }
    if(length(dbPattern) == 1) {
        dbPattern <- c(dbPattern, '\\1')
    }
    obEv <- formatObisEvent(gps(x), dbPattern = dbPattern, ...)
    es <- eventSummary(x)
    obOc <- formatObisOccur(es, speciesMap = speciesMap, dbPattern = dbPattern, ...)
    parMatch <- select(obEv, all_of(c('db', 'parentEventID' = 'eventID')))
    obOc <- left_join(obOc, parMatch, by=join_by('db'))
    obOc$db <- NULL
    obEv$db <- NULL
    list(event=obEv, occurrence=obOc)
}
