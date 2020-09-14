# random utils

# for converting from database UTC columns that are characters
pgDateToPosix <- function(x) {
    as.POSIXct(as.character(x), format='%Y-%m-%d %H:%M:%OS', tz='UTC')
}

# drop columns with names cols
dropCols <- function(x, cols) {
    keepCols <- !(colnames(x) %in% cols)
    x[, keepCols, drop=FALSE]
}

# what event is a UID in returns named index
whereUID <- function(study, UID, quiet=FALSE) {
    UID <- as.character(UID)
    where <- sapply(UID, function(u) { #for each uid
        ix <- which(sapply(events(study), function(e) { #go through each event
            any(sapply(detectors(e), function(d) { #and every detector in that event
                u %in% d$UID
            }))
        }))
        if(length(ix) == 0) {
            return(NA)
        }
        ix
    }, USE.NAMES=TRUE, simplify=FALSE)
    whereNA <- is.na(where)
    if(!quiet && any(whereNA)) {
        warning('UID(s) ', paste0(UID[whereNA], collapse=', '),
                ' not found in any events.')
    }
    where
}

# match SR function
# data needs UTC, thats it
# db is sound acq data, either DF or db
# safe to fail with NULL instead of error
#' @importFrom data.table data.table setkeyv
#'
matchSR <- function(data, db, extraCols = NULL, safe=FALSE) {
    if(is.character(db)) {
        if(!file.exists(db)) {
            if(safe) return(NULL)
            stop('Database ', db, ' does not exist.')
        }
        con <-dbConnect(db, drv=SQLite())
        on.exit(dbDisconnect(con))
        soundAcquisition <- dbReadTable(con, 'Sound_Acquisition')
        soundAcquisition$UTC <- pgDateToPosix(soundAcquisition$UTC)
    }
    if(is.data.frame(db)) {
        soundAcquisition <- db
    }
    if(!('UTC' %in% colnames(data)) ||
       !inherits(data$UTC, 'POSIXct')) {
        stop('Data must have a column "UTC" in POSIXct format.')
    }

    soundAcquisition <- soundAcquisition %>%
        mutate(Status = str_trim(.data$Status),
               SystemType = str_trim(.data$SystemType)) %>%
        filter(.data$Status=='Start') %>%
        arrange(.data$UTC) %>%
        select_(.dots = c('UTC', 'sampleRate', extraCols)) %>%
        distinct() %>%
        data.table()

    setkeyv(soundAcquisition, 'UTC')

    data <- data.table(data)
    setkeyv(data, 'UTC')

    # This rolling join rolls to the first time before. Since we filtered to only starts, it goes back
    # to whatever the last Start was.
    data <- soundAcquisition[data, roll = TRUE] %>%
        data.frame()
    srNa <- which(is.na(data$sampleRate))
    if(length(srNa) == nrow(data)) {
        srReplace <- as.integer(
            readline(prompt = 'No Sample Rate found in SoundAcquisition table. Enter Sample Rate for this data:\n')
        )
        data$sampleRate[srNa] <- srReplace
    } else if(length(srNa) > 0) {
        # get mode
        mode <- which.max(tabulate(data$sampleRate[-srNa]))
        srChoice <- menu(title=paste0('Could not get Sample Rate for all detections from the "SoundAcquistion" table.',
                                      ' Should missing values be replaced with ', mode, '(value found in table).'),
                         choices = c('Yes', 'No (I want to enter my own SR)'))
        srReplace <- switch(srChoice,
                            '1' = mode,
                            '2' = readline(prompt='What Sample Rate should be used?\n'),
                            stop('Sample Rate required for calculations.')
        )
        data$sampleRate[srNa] <- srReplace
    }
    data
}

# check if in start/stop interval
# bounds is a single start/stop, sa is sound acq table from db
#' @importFrom tidyr spread
#'
inInterval <- function(bounds, sa) {
    sa <- sa[sa$Status %in% c('Start', 'Stop'), c('UTC', 'Status', 'sampleRate')]
    first <- min(which(sa$Status == 'Start'))
    last <- max(which(sa$Status == 'Stop'))
    sa <- sa[first:last,]
    alt <- sa$Status[1:(nrow(sa)-1)] != sa$Status[2:nrow(sa)]
    sa <- sa[c(TRUE, alt), ]
    sa$id <- rep(1:(nrow(sa)/2), each=2)
    sa <- tidyr::spread(sa, 'Status', 'UTC')
    startIn <- (any((bounds[1] >= sa[['Start']]) & (bounds[1] <= sa[['Stop']])))
    endIn <- (any((bounds[2] >= sa[['Start']]) & (bounds[2] <= sa[['Stop']])))
    contain <- (any((bounds[1] <= sa[['Start']]) & (bounds[2] >= sa[['Stop']])))
    startIn || endIn || contain
}

# add list without replacing old one, only replace matching names
safeListAdd <- function(x, value) {
    if(!is.list(value) ||
       is.null(names(value))) {
        stop('Can only add named lists ')
    }
    hasName <- names(value) %in% names(x)
    if(any(hasName)) {
        for(n in names(value)[hasName]) {
            x[[n]] <- value[[n]]
        }
    }
    if(any(!hasName)) {
        x <- c(x, value[!hasName])
    }
    x
}

# named vector for AcEv, or named list of named vectors for AcSt
getEventTime <- function(x) {
    if(is.AcousticStudy(x)) {
        return(lapply(events(x), getEventTime))
    }
    allDets <- bind_rows(lapply(detectors(x), function(d) {
        d[, 'UTC', drop = FALSE]
    }))
    c(start=min(allDets$UTC), end=max(allDets$UTC))
}

# clip of fixed length, zeropads if needed and deals with edge case
clipAroundPeak <- function(wave, length) {
    if(length(wave) < length) {
        return(c(wave, rep(0, length - length(wave))))
    }
    peak <- which.max(abs(wave))
    low <- peak - floor(length/2)
    high <- ceiling(peak + length/2) -1
    if(low < 1) {
        return(wave[1:length])
    }
    if(high > length(wave)) {
        return(wave[(length(wave)-length+1):length(wave)])
    }
    wave[low:high]
}
