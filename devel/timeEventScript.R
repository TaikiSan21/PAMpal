source('timeEventFunctions.R')
library(PAMpal)
# change these to your locations
# updated 06-06 accidentally had click data
# updated 6-20 making code more flexible to different kinds of events
# updated 7-1-2022 adding possibility to work on mode='recording'
db <- '../Data/TestGPL/GPLwhat.sqlite3'
bins <- '../Data/TestGPL/Binary/'
pps <- PAMpalSettings(db=db, binaries = bins, sr_hz='auto', winLen_sec=.0025, filterfrom_khz=10, filterto_khz=NULL)
#### IF DOING BY RECORDING FILE INSTEAD SKIP TO LINE 36 NOW AND RUN mode='recording'
# this function takes time because it has to open all binaries
binTimes <- makeBinRanges(pps@binaries$list)
#### THIS MAKES EVENTS, CHANGE "6" TO WHATEVER TIME YOU NEED ####
# For doing the version of 6 minute events then splitting into 2 minute sub-events, start with 6 here
# and follow further isntructions below

# plots it makes have info on what they should look like, let me know if it looks weird
baseEvents <- makeTimeEvents(units='mins', length=6, bin=binTimes)$timeRange
baseEvents$id <- as.character(1:nrow(baseEvents))

#### RUN THIS TO SPLIT INTO 2 MIN SUBGROUPS ####
timeEv <- vector('list', length=nrow(baseEvents))
for(i in seq_along(timeEv)) {
    this <- makeTimeRanges(start=baseEvents$start[i], end=baseEvents$end[i], length=2, units='mins')
    this$id <- paste0(baseEvents$id[i], '_', 1:3)
    timeEv[[i]] <- this
}
timeEv <- bind_rows(timeEv)
timeEv$interval <- NULL
#### OTHERWISE RUN THIS ####
timeEv <- baseEvents
timeEv$id <- as.character(1:nrow(timeEv))
#### NEXT STEP ####
# process data so we can add to database
data <- processPgDetections(pps, mode='time', grouping=timeEv)
# data <- processPgDetections(pps, mode='recording')
# change 'Blue' to whatever seems appropriate for an initial label
data <- setSpecies(data, method='manual', value='Blue')
library(PAMmisc)
for(e in seq_along(events(data))) {
    thisEv <- data[[e]]
    uids <- bind_rows(lapply(getDetectorData(thisEv), function(x) {
        x[c('UTC', 'UID')]
    }))
    if(nrow(uids) == 0) {
        next
    }
    addPgEvent(db = files(thisEv)$db,
               binary = files(thisEv)$binaries,
               eventType = species(thisEv)$id,
               UIDs = unique(uids$UID),
               type = 'dg',
               comment = paste0('Added by PAMpal, event ID: ', id(thisEv)))
}

# get 20% and 10% subsets
#### CHANGE THIS TO WHATEVER TYPE OF DETECTION YOU ARE USING CLICK/GPL/ETC. ####
df <- getGPLData(data)
#-------------------------------#
#### FOR WAV RECORDING VERSION ####
detId <- unique(df$eventId)
baseEvents <- ancillary(data)$grouping
allId <- unique(baseEvents$id)
#### FOR TIME EVENT VERSION ####
detId <- unique(gsub('_[123]', '', df$eventId))
allId <- unique(gsub('_[123]', '', timeEv$id))
#--------------------------------#
# This has 20% of 6 minute files, labeled as 1_X, 2_X etc.
subset20 <- data.frame(db = gsub('\\.sqlite3', '', basename(db)),
                       evName = sample(detId, size=0.2 * length(detId), replace=FALSE))

subset20 <- left_join(subset20, baseEvents[c('start', 'end', 'id')], by=c('evName'='id'))
#### SKIP THIS _X PART IF NOT DOING SUBEVENTS FOR TIME EVENTS ####
subset20$evName <- paste0(subset20$evName, '_X')

subset20 <- arrange(subset20, start)
write.csv(subset20, file=paste0(gsub('\\.sqlite3', '', basename(db)), '_Val20.csv'))

noDetId <- allId[!allId %in% detId]
# This has 10% of files with no detections
subset10 <- data.frame(db = gsub('\\.sqlite3', '', basename(db)),
                       evName = sample(noDetId, size=0.1 * length(noDetId), replace=FALSE))
subset10 <- left_join(subset10, baseEvents[c('start', 'end', 'id')], by=c('evName'='id'))
#### SKIP THIS _X PART IF NOT DOING SUBEVENTS FOR TIME EVENTS ####
subset10$evName <- paste0(subset10$evName, '_X')

subset10 <- arrange(subset10, start)
write.csv(subset10, file=paste0(gsub('\\.sqlite3', '', basename(db)), '_NoDet10.csv'))