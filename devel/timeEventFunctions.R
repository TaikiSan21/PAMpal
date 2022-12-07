library(dplyr)
library(PamBinaries)
library(lubridate)
library(PAMpal)
library(patchwork)

# v 2.0 updated 5/31
# v 2.1 updated 6/30/22 changing the "tryFix" mode bc was bad on some data Cory had where it
# v 2.2 updated 2022-12-6 making plots with ggplot to avoid margin error
# decided to not duty cycle for a bit

makeBinRanges <- function(x, progress=TRUE) {
    tStart <- Sys.time()
    if(progress) {
        pb <- txtProgressBar(min=0, max=length(x), style=3)
        pbIx <- 0
    }
    result <- bind_rows(lapply(x, function(b) {
        b <- loadPamguardBinaryFile(b, skipLarge = TRUE, skipData = TRUE,convertDate = FALSE)
        out <- list(numStart=b$fileInfo$fileHeader$dataDate,
                    numEnd = b$fileInfo$fileFooter$dataDate)
        out$start <- convertPgDate(out$numStart)
        out$end <- convertPgDate(out$numEnd)
        out$diff <- out$numEnd - out$numStart
        out$file <- b$fileInfo$fileName
        if(progress) {
            pbIx <<- pbIx + 1
            setTxtProgressBar(pb, value=pbIx)
        }
        out
    })
    )
    result <- distinct(result, start, end, .keep_all=TRUE)
    result$interval <- interval(result$start, result$end)
    tEnd <- Sys.time()
    # cat('Time elapsed:', round(as.numeric(difftime(tEnd, tStart, units='mins')), 2), 'minutes')
    result
}

makeTimeRanges <- function(start, end=20, length='2/2', units=c('secs', 'mins', 'hours')) {
    if(is.character(length)) {
        splitLen <- strsplit(length, '/')[[1]]
        onLen <- as.numeric(splitLen[1])
        if(length(splitLen) == 1) {
            offLen <- 0
        } else {
            offLen <- as.numeric(splitLen[2])
        }
    } else {
        onLen <- length
        offLen <- 0
    }
    units <- match.arg(units)
    unitScale <- switch(units,
                        'secs' = 1,
                        'mins' = 60,
                        'hours' = 3600)
    onLen <- onLen * unitScale
    offLen <- offLen * unitScale
    if(inherits(end, 'POSIXct')) {
        totLength <- as.numeric(difftime(end, start, units='secs'))
    } else {
        totLength <- end * unitScale
    }
    startSecs <- seq(from=0, to=totLength, by=onLen + offLen)
    if(startSecs[length(startSecs)] == totLength) {
        startSecs <- startSecs[-length(startSecs)]
    }
    result <- data.frame(start = start + startSecs)
    result$end <- result$start + onLen
    result$interval <- interval(result$start, result$end)
    result
}

checkOverlap <- function(x, y, early=FALSE) {
    if(early) {
        overlap <- rep(0, nrow(x))
        for(i in seq_along(overlap)) {
            isOverlap <- int_overlaps(x$interval[i], y$interval)
            if(!any(isOverlap)) {
                overlap[i] <- 0
                next
            }
            overlap[i] <- sum(int_length(intersect(x$interval[i], y$interval[isOverlap])), na.rm=TRUE)
            if(i > 1 &&
               overlap[i] != overlap[i-1]) {
                break
            }
        }
        x$overlap <- overlap
    } else {
        x$overlap <- sapply(x$interval, function(i) {
            isOverlap <- int_overlaps(i, y$interval)
            if(!any(isOverlap)) {
                return(0)
            }
            sum(int_length(intersect(i, y$interval[isOverlap])), na.rm=TRUE)
        })
    }
    x$overlapPct <- x$overlap / int_length(x$interval)
    x
}

makeTimeEvents <- function(start=NULL, end=NULL, length, units=c('secs', 'mins', 'hours'), bin, tryFix=TRUE, plot=TRUE, progress=TRUE) {
    if(is.PAMpalSettings(bin)) {
        bin <- bin@binaries$list
    }
    if(is.character(bin)) {
        binRange <- makeBinRanges(bin, progress)
    }
    if(is.data.frame(bin)) {
        binRange <- bin
    }
    if(is.null(start)) {
        start <- min(binRange$start)
    }
    if(is.null(end)) {
        end <- max(binRange$end)
    }
    # isCont <- median(as.numeric(difftime(
    #     binRange$end, binRange$start, units='secs'
    # )))
    # isCont <- isCont < 30
    # browser()
    if(isFALSE(tryFix)) {
        timeRange <- makeTimeRanges(start, end, length=length, units=units)
    } else {
        timeRange <- fixTimeRange(start=start, end=end, bin=binRange, length=length, units=units)
    }
    binFilt <- filter(binRange,
                      end >= timeRange$start[1],
                      start <= timeRange$end[nrow(timeRange)])
    binFilt <- checkOverlap(binFilt, timeRange)
    timeRange <- checkOverlap(timeRange, binFilt)
    if(plot) {
        nOut <- nrow(binRange) - nrow(binFilt)
        n99 <- sum(binFilt$overlapPct < .99)
        tDiff <- as.numeric(difftime(timeRange$start[2:nrow(timeRange)], timeRange$end[1:(nrow(timeRange)-1)], units=units))
        # op <- par(mfrow=c(1,3))
        # on.exit(par(op))
        # hist(binFilt$overlapPct, breaks=seq(from=0, to=1, by=.02),
        #      xlim=c(0,1),
        #      main=paste0('Pct of each binary file in event ',
        #                  '\n(', nOut, ' files outside of time range, ',
        #                  n99, ' files < .99)',
        #                  '\n(All 1 unless duty cycle mismatch btwn event/recorder)'))
        g1 <- ggplot(binFilt, aes(x=overlapPct)) +
            geom_histogram(binwidth=.02) +
            xlim(-.03, 1.03) +
            labs(title=paste0('Pct of each binary file in event ',
                              '\n(', nOut, ' files outside of time range, ',
                              n99, ' files < .99)',
                              '\n(All 1 unless duty cycle mismatch btwn event/recorder)')) +
            scale_y_continuous(expand=expansion(mult=c(0, .05)))
        # hist(timeRange$overlapPct, breaks=seq(from=0, to=1, by=.02), xlim=c(0,1),
        #      main='Pct of each event with binary data\n(Should all be 1)')
        g2 <- ggplot(timeRange, aes(x=overlapPct)) +
            geom_histogram(binwidth=.02) +
            xlim(-.03, 1.03) +
            labs(title='Pct of each event with binary data\n(Should all be 1)') +
            scale_y_continuous(expand=expansion(mult=c(0, .05)))
        # if(max(tDiff) < 1) {
        #     breaks <- 'Sturges'
        # } else {
        #     breaks <- 0:ceiling(max(tDiff))
        # }
        # hist(tDiff, breaks=breaks, main=paste0('Time between events (', units, ')'), xlim=c(0, max(tDiff) + 1))
        g3 <- ggplot(data.frame(timeDiff=tDiff), aes(x=timeDiff)) +
            geom_histogram(bins=50) +
            labs(title=paste0('Time between events (', units, ')')) +
            scale_y_continuous(expand=expansion(mult=c(0, .05)))
        print(g1/g2/g3)
    }
    list(timeRange=timeRange, binRange=binFilt, allBin=binRange)
}

fixTimeRange <- function(start=NULL, end=NULL, bin, length='2/8', units='mins') {
    if(is.null(start)) {
        start <- min(bin$start)
    }
    if(is.null(end)) {
        end <- max(bin$end)
    }
    bin <- bin[(bin$end >= start) & (bin$start <= end), ]
    time <- makeTimeRanges(start=start, end=end, length=length, units=units)
    time <- checkOverlap(time, bin, early=TRUE)
    # isCont <- median(as.numeric(difftime(
    #     bin$start, bin$end, units='secs'
    # )))
    # isCont <- isCont < 30
    # # if recs are cont we dont need to reset all the time
    # # it will sort itself out
    # if(isCont) {
    #     return(time)
    # }
    result <- list()
    newBin <- bin
    newTime <- time
    pb <- txtProgressBar(min=0, max=nrow(newTime), style=3)
    for(i in 1:nrow(time)) {
        # cat(paste0(nrow(newTime), ' rows remaining...\n'))
        if(nrow(newTime) <= 1) {
            result[[i]] <- newTime
            setTxtProgressBar(pb, value=nrow(time))
            break
        }
        # browser()
        change <- which(newTime$overlap[2:nrow(newTime)] != newTime$overlap[1:(nrow(newTime)-1)])
        if(length(change) == 0) {
            result[[i]] <- newTime
            setTxtProgressBar(pb, value=nrow(time))
            break
        }
        change <- min(change) #+ 1
        result[[i]] <- newTime[1:change, ]
        lastEnd <- newTime$end[change+1]
        lastStart <- newTime$start[change+1]
        endIn <- lastEnd %within% newBin$interval &
            lastEnd != newBin$start
        startIn <- lastStart %within% newBin$interval &
            lastStart != newBin$end
        if(any(startIn)) {
            startIx <- min(which(startIn))
        } else if(any(endIn)) {
            startIx <- min(which(endIn))
        } else {
            startIx <- min(which(newBin$start > lastEnd))
        }
        newBin <- newBin[startIx:nrow(newBin), ]
        if(nrow(newBin) == 0) {
            break
        }
        nextStart <- min(newBin$start)
        # if(any(endIn)) {
        #     newBin <- newBin[min(which(endIn)):nrow(newBin), ]
        #     if(nrow(newBin) == 0) {
        #         break
        #     }
        #     nextStart <- min(newBin$start)
        # } else if(any(startIn)) {
        #     newBin <- newBin[min(which(startIn)):nrow(newBin), ]
        #     if(nrow(newBin) == 0) {
        #         break
        #     }
        #     nextStart <- min(newBin$start)
        # } else {
        #     newBin <- newBin[newBin$start > lastEnd, ]
        #     if(nrow(newBin) == 0) {
        #         break
        #     }
        #     nextStart <- min(newBin$start)
        # }
        nextStart <- max(nextStart, newTime$end[change])
        # nextBinTime <- min(newBin$start)
        newTime <- makeTimeRanges(start=nextStart, end=max(newBin$end), length=length, units=units)
        newTime <- checkOverlap(newTime, newBin, early=TRUE)
        if(i/nrow(time) > .79) {
            # browser()
        }
        setTxtProgressBar(pb, value = nrow(time) - nrow(newTime))
    }
    bind_rows(result)
}