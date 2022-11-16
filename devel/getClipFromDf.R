getClipFromDf <- function(df, rec, buffer = c(0, 0.1), mode=c('detection'),
                          channel = 1, useSample=FALSE, progress=TRUE, verbose=TRUE, FUN=NULL, ...) {
    if(is.character(rec)) {
        wavMap <- PAMpal:::mapWavFolder(rec)
    } else {
        wavMap <- rec
    }
    wavMap$numStart <- as.numeric(wavMap$start)
    wavMap$numEnd <- as.numeric(wavMap$end)
    # mode <- match.arg(mode)


    # dbMap <- split(recs, recs$db)
    # names(dbMap) <- basename(names(dbMap))
    #
    if(length(buffer) == 1) {
        buffer <- buffer * c(-1, 1)
    }
    buffer <- abs(buffer) * c(-1, 1)

    if(is.null(FUN)) {
        FUN <- function(wav, ...) wav
    }
    # each database can have a different set of assigned recordings,
    # so we break up by DB

    # setup warning storage
    noMatch <- character(0)
    multMatch <- character(0)
    nonConsec <- character(0)
    fileDNE <- character(0)
    noChan <- character(0)
    result <- character(0)
    on.exit({
        if(length(noMatch) > 0) {
            warning('Could not find matching wav files for ', mode, ' ', PAMpal:::printN(noMatch, 6), call.=FALSE)
        }
        if(length(multMatch) > 0) {
            warning(PAMpal:::oneUpper(mode), ' ', PAMpal:::printN(multMatch, 6),
                    ' matched more than 1 possible wav file, could not get clip.', call.=FALSE)
        }
        if(length(nonConsec) > 0) {
            warning(PAMpal:::oneUpper(mode), ' ', PAMpal:::printN(nonConsec, 6),
                    ' spanned two non-consecutive wav files, could not get clip.', call.=FALSE)
        }
        if(length(fileDNE) > 0) {
            warning('Wav files for ', mode, ' ', PAMpal:::printN(fileDNE, 6), ' could not be found on disk.',
                    ' Function "updateFiles" can help relocate files that have moved.', call. = FALSE)
        }
        if(length(noChan) > 0) {
            warning('Wav files for ', mode, ' ', PAMpal:::printN(noChan, 6),
                    ' did not have the desired channels.', call.=FALSE)
        }
    })
    # # one DB at a time
    # for(d in seq_along(dbMap)) {
    # thisDbData <- x[which(evDbs == names(dbMap)[d])]
    # if(length(events(thisDbData)) == 0) next
    # allTimes <- getTimeRange(thisDbData, mode=mode, sample=useSample)
    allResult <- vector('list', length = nrow(df))
    names(allResult) <- df$id
    if(progress) {
        cat('Processing wav files...\n', sep='')
        pb <- txtProgressBar(min=0, max=length(allResult), style = 3)
    }
    for(i in seq_along(allResult)) {
        # timeRange <- c(allTimes[[i]]$start, allTimes[[i]]$end) + buffer
        timeRange <- c(df$start[i], df$end[i]) + buffer
        # for start and end check if in range. if we buffered, try undoing that first.
        # so like if buffer put us before first file, start and beginning of first file instead.
        startIx <- PAMpal:::checkIn(timeRange[1], wavMap)
        if(length(startIx) > 1) {
            multMatch <- c(multMatch, names(allResult)[i])
            if(progress) {
                setTxtProgressBar(pb, value=i)
            }
            next
        }
        if(is.na(startIx)) {
            startIx <- PAMpal:::checkIn(timeRange[1] - buffer[1], wavMap)
            if(is.na(startIx)) {
                noMatch <- c(noMatch, names(allResult)[i])
                if(progress) {
                    setTxtProgressBar(pb, value=i)
                }
                next
            }
            timeRange[1] <- wavMap$start[startIx]
        }
        endIx <- PAMpal:::checkIn(timeRange[2], wavMap)
        if(length(endIx) > 1) {
            multMatch <- c(multMatch, names(allResult)[i])
            if(progress) {
                setTxtProgressBar(pb, value=i)
            }
            next
        }
        if(is.na(endIx)) {
            endIx <- PAMpal:::checkIn(timeRange[2] - buffer[2], wavMap)
            if(is.na(endIx)) {
                noMatch <- c(noMatch, names(allResult)[i])
                if(progress) {
                    setTxtProgressBar(pb, value=i)
                }
                next
            }
            timeRange[2] <- wavMap$end[endIx]
        }
        if(wavMap$fileGroup[startIx] != wavMap$fileGroup[endIx]) {
            nonConsec <- c(nonConsec, names(allResult)[i])
            if(progress) {
                setTxtProgressBar(pb, value=i)
            }
            next
        }
        if(any(!file.exists(wavMap$file[startIx:endIx]))) {
            fileDNE <- c(fileDNE, names(allResult)[i])
            if(progress) {
                setTxtProgressBar(pb, value=i)
            }
            next
        }
        startTime <- as.numeric(difftime(timeRange[1], wavMap$start[startIx], units='secs'))
        endTime <- as.numeric(difftime(timeRange[2], wavMap$start[endIx], units='secs'))
        wavResult <- vector('list', length = endIx)
        for(w in startIx:endIx) {
            readStart <- 0
            readEnd <- Inf
            if(w == startIx) {
                readStart <- startTime
            }
            if(w == endIx) {
                readEnd <- endTime
            }
            wavResult[[w]] <- tuneR::readWave(wavMap$file[w], from = readStart, to = readEnd, units = 'seconds', toWaveMC = TRUE)
        }

        wavResult <- wavResult[!sapply(wavResult, is.null)]
        wavResult <- do.call(bind, wavResult) # [, 1:min(2, ncol(wavResult))]
        chanIn <- channel <= ncol(wavResult@.Data)
        if(!any(chanIn)) {
            noChan <- c(noChan, names(allResult)[i])
            if(progress) {
                setTxtProgressBar(pb, value = i)
            }
            next
        }
        if(!all(chanIn)) {
            noChan <- c(noChan, names(allResult)[i])
        }
        wavResult <- wavResult[, channel[chanIn]]
        colnames(wavResult) <- MCnames$name[1:ncol(wavResult)]
        allResult[[i]] <- FUN(wavResult, name=names(allResult)[i], time=timeRange, channel=channel[chanIn], mode=mode, ...)

        if(progress) {
            setTxtProgressBar(pb, value=i)
        }
    }
    isNull <- sapply(allResult, is.null)
    if(verbose) {
        cat('\n', paste0('Processed ', sum(!isNull), ' wav file(s).\n'))
    }
    result <- c(result, allResult)
    # }
    invisible(result)
}
