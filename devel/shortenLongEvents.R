library(lubridate)
# updated 2023-01-11 runs faster
longEvents <- checkStudy(myData_pascal, maxLength = 120)#168events

getTimes <- function(event) {
    # get all the detector data in this event
    allDets <- getDetectorData(event)
    # This goes through all of the $click, $whistle, and $cepstrum detectors if present
    justTimes <- bind_rows(lapply(allDets, function(x) {
        x[, c('UID', 'UTC')]
    }))
    justTimes
}

timeToStartEnd <- function(time, length = 120) {
    range <- range(time)
    lenSecs <- as.numeric(difftime(range[2], range[1], units='secs'))
    # Figure out how many events of length "length" we can have
    numSplits <- ceiling(lenSecs / length)
    # Each event starts some multiple of length from the original start
    start <- range[1] + length * 0:(numSplits-1)
    # Doesnt matter that the end is different than the actual event end
    # here since we are just using these times to filter later
    end <- start + length
    list(start=start, end=end, interval=interval(start, end))
}

withinLHS <- function(a, int) {
    int <- int_standardize(int)
    as.numeric(a) - as.numeric(int@start) < int@.Data & as.numeric(a) - as.numeric(int@start) >= 0
}

newEvents <- vector('list', length = length(events(myData_pascal)))
pb <- txtProgressBar(min=0, max=length(newEvents), style=3)
for(i in seq_along(newEvents)) {
    thisEvent <- myData_pascal[[i]]
    # leftovers from old versions of PAMpal take up space
    if('source' %in% names(settings(thisEvent))) {
        settings(thisEvent)$source <- NULL
    }
    # removing ICI because these should be recalculated after split
    if('ici' %in% names(ancillary(thisEvent))) {
        thisEvent@ancillary$ici <- NULL
    }
    thisTime <- getTimes(thisEvent)
    thisStartEnd <- timeToStartEnd(thisTime$UTC)
    if(length(thisStartEnd$start) == 1) {
        newEvents[[i]] <- list(thisEvent)
        setTxtProgressBar(pb, value=i)
        next
    } 
    evList <- vector('list', length = length(thisStartEnd$start))
    for(s in seq_along(thisStartEnd$start)) {
        onePart <- thisEvent
        dets <- lapply(detectors(onePart), function(x) {
            x <- x[withinLHS(x$UTC, thisStartEnd$interval[s]), ]
            if(nrow(x) == 0) {
                return(NULL)
            }
            x
        })
        dets <- dets[sapply(dets, function(x) !is.null(x))]
        if(length(dets) == 0) next

        detectors(onePart) <- dets
        id(onePart) <- paste0(id(onePart), '_', s)
        evList[[s]] <- onePart
    }
    newEvents[[i]] <- evList
    setTxtProgressBar(pb, value=i)
}

newEvents <- unlist(newEvents)
names(newEvents) <- sapply(newEvents, id)
shortData <- myData_pascal
events(shortData) <- newEvents
myData_pascal #see how many events are here
shortData #see how many there are now
noWarns <- checkStudy(shortData, maxLength = 120)
noWarns
