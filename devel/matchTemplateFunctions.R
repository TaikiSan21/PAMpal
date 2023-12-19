library(dplyr)
library(PamBinaries)
library(PAMmisc)

loadTemplateBinary <- function(x, names, columns) {
    bin <- loadPamguardBinaryFile(x, skipLarge=TRUE, convertDate=FALSE)
    if(length(bin$data) == 0) {
        return(NULL)
    }
    bin <- try(suppressWarnings(pbToDf(bin, templateNames=names)[columns]))
    if(inherits(bin, 'try-error')) {
        stop('Problem loading binary file ', x)
    }
    bin$UID <- as.character(bin$UID)
    bin$BinaryFile <- basename(x)
    bin$date <- convertPgDate(bin$date)
    if(!is.data.frame(bin)) {
        warning('Weird bin')
        print(bin)
        return(NULL)
    }
    bin <- rename(bin, UTC = date)
    bin
}

loadTemplateFolder <- function(dir, names, extraCols=NULL, file=NULL, progress=TRUE) {
    if(!dir.exists(dir)) {
        stop('Directory ', dir, ' does not exists')
    }
    files <- list.files(dir, pattern='^Click_Detector.*pgdf$', full.names=TRUE, recursive=TRUE)
    columns <- c('UID', 'type', 'date',
                 paste0(names, '_thresh'),
                 extraCols)
    if(progress) {
        pb <- txtProgressBar(min=0, max=length(files), style=3)
        ix <- 0
    }
    allBins <- bind_rows(lapply(files, function(x) {
        if(progress) {
            ix <<- ix + 1
            setTxtProgressBar(pb, value=ix)
        }
        loadTemplateBinary(x, names=names, columns=columns)
    }))
    noTemplate <- is.na(allBins[[paste0(names[1], '_match')]])
    if(any(noTemplate)) {
        warning(sum(noTemplate), ' detections did not have template classifier data')
    }
    if(!is.null(file)) {
        saveRDS(allBins, file=file)
    }
    allBins
}

addTemplateLabels <- function(x, db, names, thresh) {
    if(!file.exists(db)) {
        stop('Database ', db, ' does not exist')
    }
    eventData <- PAMpal:::getDbData(db, doSR=FALSE)
    labelCols <- c('parentID', 'eventLabel', 'templateId', 'templateMatch')
    x <- PAMpal:::dropCols(x, labelCols)
    if(is.null(eventData) || nrow(eventData) == 0) {
        x$parentID <- NA
        x$eventLabel <- NA
    } else {
        eventData <- distinct(eventData[c('UID', 'parentID', 'eventLabel', 'BinaryFile')])
        x <- left_join(x, eventData,
                       by=join_by('UID'=='UID', 'BinaryFile'=='BinaryFile'))
    }
    noMatch <- is.na(x$parentID)
    x$parentID[noMatch] <- 'none'
    x$eventLabel[noMatch] <- 'none'
    threshCols <- paste0(names, '_thresh')
    if(length(thresh) != length(threshCols)) {
        stop('Number of threshold values must match number of template names')
    }
    aboveThresh <- rep('', nrow(x))
    for(i in seq_along(threshCols)) {
        isAbove <- x[[threshCols[i]]] > thresh[i]
        isAbove[is.na(isAbove)] <- FALSE
        aboveThresh <- paste0(aboveThresh,
                              ifelse(isAbove, i, ''))
    }
    aboveThresh[aboveThresh == ''] <- 'none'
    x$templateId <- aboveThresh
    x$templateMatch <- aboveThresh != 'none'
    x
}

markGoodEvents <- function(x, minDets=3, maxSep=120, maxLength=120, nDets = 3, nSeconds=120, maxSeconds=120, verbose=TRUE) {
    x <- PAMpal:::dropCols(x, c('templateEvent'))
    if(!missing(nDets)) {
        warning('"nDets" has been renamed to "minDets" for clarity')
        minDets <- nDets
    }
    if(!missing(nSeconds)) {
        warning('"nSeconds" has been renamed to "maxSep" for clarity')
        maxSep <- nSeconds
    }
    if(!missing(maxSeconds)) {
        warning('"maxSeconds" has been renamed to "maxLength" for clarity')
        maxLength <- maxSeconds
    }
    result <- bind_rows(
        lapply(
            split(x, x$templateMatch), function(d) {
                if(isFALSE(d$templateMatch[1])) {
                    d$templateEvent <- 'none'
                    return(d)
                }
                d <- arrange(d, UTC)
                d$tDiff <- c(0, as.numeric(difftime(d$UTC[2:nrow(d)],
                                                    d$UTC[1:(nrow(d)-1)], units='secs')))
                d$overtime <- d$tDiff > maxSep
                d$templateEvent <- as.character(cumsum(d$overtime))
                eventDets <- table(d$templateEvent)
                goodEvents <- names(eventDets)[eventDets >= minDets]
                d$templateEvent[!d$templateEvent %in% goodEvents] <- 'none'
                d <- PAMpal:::dropCols(d, c('tDiff', 'overtime'))
                d
            }))
    # SHOULD REALLY RESET START TIMER AFTER EVERY NEW EVENT
    result <- bind_rows(
        lapply(
            split(result, result$templateEvent), function(e) {
                if(e$templateEvent[1] == 'none') {
                    return(e)
                }
                tDiff <- c(0, as.numeric(difftime(e$UTC[2:nrow(e)],
                                                    e$UTC[1:(nrow(e)-1)], units='secs')))
                # e$tInterval <- floor(cumsum(e$tDiff) / maxLength)
                thisLab <- 1
                labs <- rep(thisLab, nrow(e))
                cTimes <- cumsum(tDiff)
                for(i in 1:nrow(e)) {
                    if(cTimes[i] <= maxLength) {
                        next
                    }
                    thisLab <- thisLab + 1
                    labs[i:(length(labs))] <- thisLab
                    cTimes <- cTimes - cTimes[i]
                }
                e$templateEvent <- paste0(e$templateEvent, '_', labs)
                e
            }
        ))
    result$truePos <- (result$parentID != 'none') & (result$templateEvent != 'none')
    result$trueNeg <- (result$parentID == 'none') & (result$templateEvent == 'none')
    result$falsePos <- (result$parentID == 'none') & (result$templateEvent != 'none')
    result$falseNeg <- (result$parentID != 'none') & (result$templateEvent == 'none')
    if(verbose) {
        cat('TP: ', sum(result$truePos), '\n',
            'FP: ', sum(result$falsePos), '\n',
            'TN: ', sum(result$trueNeg), '\n',
            'FN: ', sum(result$falseNeg), '\n', sep='')
    }
    result
}

addTemplateEvents <- function(db, binFolder, data) {
    if(!file.exists(db)) {
        stop('Database ', db, ' does not exist')
    }
    if(!dir.exists(binFolder)) {
        stop('Binary folder ', binFolder, ' does not exist')
    }
    binFiles <- list.files(binFolder, pattern='^Click_Detector.*pgdf$', full.names=TRUE, recursive=TRUE)
    data <- data[data$templateEvent != 'none', ]
    data <- split(data, data$templateEvent)

    pb <- txtProgressBar(min=0, max=length(data), style=3)
    for(e in seq_along(data)) {
        bins <- binFiles[basename(binFiles) %in% unique(data[[e]]$BinaryFile)]
        label <- ifelse(any(data[[e]]$truePos), 'TP', 'FP')
        addPgEvent(db=db, binary=bins, UIDs = unique(data[[e]]$UID), type='click',
                   eventType = label,
                   comment=paste0('templateEvent: ', data[[e]]$templateEvent[1],
                                  ', posRate: ', round(mean(data[[e]]$truePos), 2)))
        setTxtProgressBar(pb, value=e)
    }
    TRUE
}

summariseManualEvents <- function(x) {
    result <- x %>%
        filter(parentID != 'none') %>%
        group_by(parentID) %>%
        summarise(eventLabel = unique(eventLabel),
                  nTemplate = sum(templateEvent != 'none'),
                  pctTemplate = mean(templateEvent != 'none'),
                  nDets = n())
    nZero <- sum(result$nTemplate == 0)
    cat('\n', nZero, ' out of ', nrow(result), ' manual event(s) had no template detections (FN).', sep='')
    result
}

summariseTemplateEvents <- function(x) {
    result <- x %>%
        filter(templateEvent != 'none') %>%
        group_by(templateEvent) %>%
        summarise(manualEvent = paste0(unique(parentID), collapse=', '),
                  manualLabel = paste0(unique(eventLabel), collapse=', '),
                  nManual = sum(parentID != 'none'),
                  pctManual = mean(parentID != 'none'),
                  nDets = n()) %>%
        # filter(oldEv == 'none') %>%
        mutate(manualEvent = ifelse(manualEvent == 'none', 'none', gsub(', none|none, ', '', manualEvent)),
               manualLabel = ifelse(manualLabel == 'none', 'none', gsub(', none|none, ', '', manualLabel)))
    nZero <- sum(result$nManual == 0)
    cat('\n', nZero, ' out of ', nrow(result), ' tempalte event(s) had no manual detections (FP).', sep='')
    multiMatch <- grepl(',', result$manualEvent)
    if(any(multiMatch)) {
        cat('\n', sum(multiMatch), ' template event(s) matched multiple manual events.', sep='')
    }
    result
}
