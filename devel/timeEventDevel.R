library(dplyr)
library(PamBinaries)
library(lubridate)
library(PAMpal)

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

checkOverlap <- function(x, y) {
    x$overlap <- sapply(x$interval, function(i) {
        isOverlap <- int_overlaps(i, y$interval)
        if(!any(isOverlap)) {
            return(0)
        }
        sum(int_length(intersect(i, y$interval[isOverlap])), na.rm=TRUE)
    })
    x$overlapPct <- x$overlap / int_length(x$interval)
    x
}

makeTimeEvents <- function(start=NULL, end=NULL, length, units=c('secs', 'mins', 'hours'), bin, plot=TRUE, progress=TRUE) {
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
    timeRange <- makeTimeRanges(start, end, length=length, units=units)
    binFilt <- filter(binRange,
                       end >= timeRange$start[1],
                       start <= timeRange$end[nrow(timeRange)])
    binFilt <- checkOverlap(binFilt, timeRange)
    timeRange <- checkOverlap(timeRange, binFilt)
    nOut <- nrow(binRange) - nrow(binFilt)
    if(plot) {
        op <- par(mfrow=c(2,1))
        on.exit(par(op))
        hist(binFilt$overlapPct, breaks=seq(from=0, to=1, by=.02),
             xlim=c(0,1), main=paste0('Percent coverage of binary files \n(', nOut, ' files outside of time range)'))
        hist(timeRange$overlapPct, breaks=seq(from=0, to=1, by=.02), xlim=c(0,1), main='Percent coverage of events')
    }
    list(timeRange=timeRange, binRange=binFilt, allBin=binRange)
}


int_length(intersect(timeRange$interval[1], binTimes$interval[1]))


fixedTimeEvents <- function(pps, startTime=NULL, length='2/18') {
}
hm <- tev$timeRange
hm$ovChange <- FALSE
hm$ovChange[2:nrow(hm)] <- hm$overlap[1:(nrow(hm)-1)] != hm$overlap[2:nrow(hm)]
sa <- PAMpal:::readSa(pps@db)
dbTimes <- sa[sa$Status == 'Start',]
binTimes <- bind_rows(lapply(pps@binaries$list, binRange))
plot(x=dbTimes, y=rep(1, length(dbTimes)))
points(x=binTimes$start, y=rep(1.2, nrow(binTimes)))

binTimes$toNext <- NA
binTimes$toNext[1:(nrow(binTimes)-1)] <- binTimes$numStart[2:nrow(binTimes)] - binTimes$numEnd[1:(nrow(binTimes)-1)]

# problems - even if continuous, pamguard has hiccups. So if you think you have 2/2 continuous, you will
# eventually glitch and your 2/2 will not line up to your data. At first point of glitch and at
# end will have poorly cut events.
