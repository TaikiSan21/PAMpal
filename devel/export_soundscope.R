#' @title Export a Soundscope NetCDF
#'
#' @description Exports the data in an \linkS4class{AcousticStudy}
#'   object to a NetCDF file compatible with the Soundscope program
#'
#' @param x a \linkS4class{AcousticStudy} object
#'
#' @return file path to the exported NetCDF file
#'
#' @author Taiki Sakai \email{taiki.sakai@@noaa.gov}
#'
#' @importFrom uuid UUIDgenerate
#' @export
#'
##################################################################
##### !!!STOP!!! UPDATE THIS FROM SOUNDSCOPEFUNCTIONS.R FIRST ####
##################################################################
export_soundscope <- function(x, detector=c('click', 'whistle', 'gpl', 'cepstrum'),
                              bin=NULL, extraCols=NULL,
                              offset=0, channel=NULL, filename) {
    recData <- files(x)$recordings
    if(is.null(recData)) {
        stop('Recording information not found, use "addRecordings" first.')
    }
    detector <- match.arg(detector)
    x <- getDetectorData(x, measures=TRUE)[[detector]]
    if('Channel' %in% colnames(x)) {
        if(is.null(channel)) {
            channel <- x$Channel[1]
        }
        x <- filter(x, Channel == channel)
    }
    if(!is.null(bin)) {
        x$UTC <- floor_date(x$UTC, unit=bin)
        x <- x %>%
            group_by(.data$UTC, .data$Channel) %>%
            summarise(across(any_of(extraCols), median),
                      n=n(),
                      species=unique(.data$species)) %>%
            ungroup()
        extraCols <- unique(c('n', extraCols))
    }

    x <- joinRecordingData(x, recData, columns=c('file', 'start', 'sr'))
    noWavMatch <- is.na(x$file)
    if(all(noWavMatch)) {
        stop('No detections had matching wav files, could not create NetCDF')
    }
    if(any(noWavMatch)) {
        warning(sum(noWavMatch), ' detections did not have matching wav files,',
                ' they will not be included in NetCDF output.')
        x <- x[!noWavMatch, ]
    }

    x <- arrange(x, 'UTC')
    nDets <- nrow(x)
    if(length(unique(x$UTC)) != nDets) {
        # stop('Detections must have unique millisecond time for NetCDF')
    }
    ncAttributes <- list()
    ncData <- list()
    ncAttributes$measurer_name <- 'PAMGuard'
    ncAttributes$measurer_version <- '1.0'
    # do these just let it look at them?
    ncAttributes$measurements_name <- c('frequency_bandwidth', 'amplitude')
    ncData$uuid <- uuid::UUIDgenerate(n=nDets)
    # just in case theres a dupe just try again
    if(length(unique(ncData$uuid)) < nDets) {
        ncData$uuid <- uuid::UUIDgenerate(n=nDets)
    }
    duration <- ifelse(is.null(bin), .0025, as.numeric(unitToPeriod(bin)))
    ncData$time_min_offset <- as.numeric(x$UTC) - as.numeric(x$start) # "startSeconds" so should be time into file?
    ncData$time_max_offset <- ncData$time_min_offset + duration # + duration
    ncData$audio_file_dir <- dirname(x$file) # join it earlier
    ncData$audio_file_name <- tools::file_path_sans_ext(basename(x$file)) # join it earlier
    ncData$audio_file_extension <- paste0('.', tools::file_ext(basename(x$file))) # join it earlier
    ncData$audio_file_start_date <- x$start
    ncData$time_min_date <- x$UTC # this was file start + startSeconds but why not just UTC?
    ncData$time_max_date <- x$UTC + duration # + dur again
    ncData$label_class <- x$species #this was manual + SpeciesCode
    ncData$duration <- duration #max - min offset. why is this not just duration???
    ncData$from_detector <- TRUE
    ncData$audio_sampling_frequency <- x$sr # decimator or wav file sr ?
    ncData$operator_name <- 'unknown'
    ncData$UTC_offset <- offset
    ncData$confidence <- x$n #1
    ncData$software_name <- detector #just do detector?
    ncData$audio_channel <- as.numeric(x$Channel) #
    ncData$entry_date <- '' # PCTime column from database. Hm.
    ncData$frequency_min <- 200 #(x$centerkHz_10dB - x$BW_10dB/2)*1e3#'' # these are just detector settings im pretty sure
    ncData$frequency_max <- 16e3 #(x$centerkHz_10dB + x$BW_10dB/2)*1e3
    ncData$frequency_bandwidth <- ''
    ncData$amplitude <- '' # never knew where this came from
    ncAttributes$deployment_file <- '' # deployment_info.csv ???
    ## adjusting UTC offset and changing utc offset???
    ncData$date <- ncData$time_min_date + offset*3600
    # these were missing ####
    ncData$software_version <- as.character(packageVersion('PAMpal'))
    ncData$audio_bit_depth <- 16
    ncData$mooring_platform_name <- ''
    ncData$recorder_type <- ''
    ncData$recorder_SN <- ''
    ncData$hydrophone_model <- ''
    ncData$hydrophone_SN <- ''
    ncData$hydrophone_depth <- 10
    ncData$location_name <- ''
    ncData$location_lat <- 40
    ncData$location_lon <- -120
    ncData$location_water_depth <- 100
    ncData$deployment_ID <- ''
    ncData$label_subclass <- ''
    if(!is.null(extraCols)) {
        for(i in seq_along(extraCols)) {
            if(extraCols[i] %in% colnames(x)) {
                ncData[[extraCols[i]]] <- x[[extraCols[i]]]
            }
        }
    }
    # ####
    ncData <- bind_rows(ncData)
    list(data=ncData, attributes=ncAttributes)
}

joinRecordingData <- function(x, recData, columns=c('file', 'start', 'sr')) {
    if(!'UTC' %in% colnames(x)) {
        stop('"x" must have column "UTC"')
    }
    matchIx <- sapply(x$UTC, function(t) {
        PAMpal:::checkIn(t, recData)[1]
    })
    if(all(is.na(matchIx))) {
        for(c in columns) {
            x[[c]] <- NA
        }
        return(x)
    }
    recData <- select(recData, any_of(columns))
    cbind(x, recData[matchIx, ])
}

createTimeVar <- function(x) {
    startTime <- x[1]
    startNum <- as.numeric(startTime)
    startMillis <- round(startNum - floor(startNum), 3) * 1e3
    x <- as.numeric(x) - startNum
    x <- round(x * 1e3, 0) # now milliseconds since start
    startString <- format(
        as.POSIXct(floor(startNum), origin='1970-01-01 00:00:00', tz='UTC'),
        format='%Y-%m-%d %H:%M:%S'
    )
    startString <- paste0(startString, '.', startMillis)
    units <- paste0('milliseconds since ', startString)
    list(vals=x, units=units)
}

createSoundscopeNc <- function(data, attributes, file='test.nc') {
    nc <- create.nc(file, prefill=FALSE, format='netcdf4', diskless=FALSE)
    dateVar <- createTimeVar(data$date)
    dim.def.nc(nc, 'date', dimlength=length(dateVar$vals))
    var.def.nc(nc, varname='date', vartype='NC_INT64', dimensions='date')
    att.put.nc(nc, 'date', 'units', 'NC_CHAR', dateVar$units)
    att.put.nc(nc, 'date', 'calendar', 'NC_CHAR', 'proleptic_gregorian')
    var.put.nc(nc, variable='date', data=dateVar$vals)
    on.exit(close.nc(nc))
    data$date <- NULL
    for(col in colnames(data)) {
        # cat('Column ', col, '\n')
        thisClass <- class(data[[col]])[1]
        if(col %in% c('audio_channel', 'audio_sampling_frequency', 'audio_bit_depth')) {
            thisClass <- 'integer'
        }
        switch(thisClass,
               'character' = {
                   var.def.nc(nc, varname=col, vartype='NC_STRING', dimensions='date')
                   var.put.nc(nc, variable=col, data=data[[col]])
               },
               'POSIXct' = {
                   dateVar <- createTimeVar(data[[col]])
                   var.def.nc(nc, varname=col, vartype='NC_INT64', dimensions='date')
                   att.put.nc(nc, variable=col, 'units', 'NC_CHAR', dateVar$units)
                   att.put.nc(nc, variable=col, 'calendar', 'NC_CHAR', 'proleptic_gregorian')
                   var.put.nc(nc, variable=col, data=dateVar$vals)
               },
               'numeric' = {
                   var.def.nc(nc, varname=col, vartype='NC_DOUBLE', dimensions='date')
                   var.put.nc(nc, variable=col, data=data[[col]])
               },
               'integer' = {
                   var.def.nc(nc, varname=col, vartype='NC_INT', dimensions='date')
                   var.put.nc(nc, variable=col, data=data[[col]])
               },
               'logical' = {
                   var.def.nc(nc, varname=col, vartype='NC_BYTE', dimensions='date')
                   att.put.nc(nc, variable=col, 'dtype', 'NC_CHAR', 'bool')
                   var.put.nc(nc, variable=col, data=as.numeric(data[[col]]))
               },
               {
                   warning('IDK what class ', thisClass, ' is')
               }
        )
    }
    attributes$datatype <- 'Measurement'
    for(att in names(attributes)) {
        att.put.nc(nc, variable='NC_GLOBAL', name=att, type='NC_CHAR', value=attributes[[att]])
    }
    print.nc(nc)
    file
}

unitToPeriod <- function(x) {
    x <- gsub('([0-9]*)(.*)', '\\1_\\2', x)
    x <- strsplit(x, '_')[[1]]
    if(x[1] == '') {
        x[1] <- '1'
    }
    period(as.numeric(x[1]), units=x[2])
}
