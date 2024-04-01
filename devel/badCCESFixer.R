cces <- readRDS('../Data/BadCCES/ccesStudy.rds')
hm <- updateFiles(cces, db='../Data/BadCCES/CCES Beaked Whale Databases w-GPS/')

hm <- addGps(hm)
View(getWarnings(hm))


plotGps <- gps(hm)
ggplot(plotGps, aes(x=Longitude, y=Latitude, col=basename(db))) +
    geom_point()
#EventInfo
load('../Data/BadCCES/CCES2018_BW_Detections.rda', verbose=T)
str(EventInfo)
# hm we have species differnces like ?bw does not exist in my CCES
#CCES2018
load('../Data/BadCCES/CCES2018_Deployment_Metadata.rda', verbose=T)
str(CCES2018)
CCES2018 <- CCES2018[!is.na(CCES2018$Drift.), ]
CCES2018$Recovery_Longitude[3] <- '128.3174'
CCES2018$Recovery_Latitude[3] <- '34.3774'
CCES2018$Recovery_Latitude[9] <- '28.2857'
CCES2018$Deployment_Longitude[3] <- '126.6449'
CCES2018$Deployment_Longitude[4] <- '125.0584'
CCES2018$Recovery_Longitude <- as.numeric(CCES2018$Recovery_Longitude)
CCES2018$Deployment_Longitude <- as.numeric(CCES2018$Deployment_Longitude)
CCES2018$Recovery_Latitude <- as.numeric(CCES2018$Recovery_Latitude)
CCES2018$Deployment_Longitude <- CCES2018$Deployment_Longitude * -1
CCES2018$Recovery_Longitude <- CCES2018$Recovery_Longitude * -1
#------CCES2018_Tracks
load('../Data/BadCCES/CCES2018_DriftTracks.rda', verbose=T)
str(CCES2018_Tracks)
table(CCES2018_Tracks$DeploymentID)
table(as.numeric(gsub('.*Drift-([0-9]{1,2}).*', '\\1', basename(plotGps$db))))
library(ggplot2)
plotGps$DeploymentID <- as.numeric(gsub('.*Drift-([0-9]{1,2}).*', '\\1', basename(plotGps$db)))
plotGps$UTCMilliseconds <- NULL
plotGps$Project <- 'CCES'
plotGps$db <- NULL
data.table::setDF(plotGps)
#------GOOD
str(plotGps)


# TESTING for track plots
library(patchwork)

outs <- vector('list', length=length(unique(CCES2018_Tracks$DeploymentID)))
names(outs) <- as.character(unique(CCES2018_Tracks$DeploymentID))
for(i in unique(CCES2018_Tracks$DeploymentID)) {
    old <- ggplot(filter(CCES2018_Tracks, DeploymentID == i),
                  aes(x=Longitude, y=Latitude)) +
        geom_point() +
        ggtitle(paste0('Old ', i))
    new <- ggplot(filter(plotGps, DeploymentID == i),
                  aes(x=Longitude, y=Latitude)) +
        geom_point()
    outs[[as.character(i)]] <- old + new
}
outs
# TESTING for CCES2018 - looks like from depDetails

for(i in unique(CCES2018_Tracks$DeploymentID)) {
    plotters <- ggplot() +
        geom_point(data=filter(plotGps, DeploymentID == i),
                   aes(x=Longitude, y=Latitude)) +
        geom_point(data=filter(CCES2018, DeploymentID == i),
                   aes(x=Deployment_Longitude, y=Deployment_Latitude), col='blue') +
        geom_point(data=filter(CCES2018, DeploymentID == i),
                   aes(x=Recovery_Longitude, y=Recovery_Latitude), col='red')
    outs[[as.character(i)]] <- plotters
}
outs
# track 8 has maybe wrog recovery location

#TESTING for EventInfo
source('devel/OBIS_Functions.R')
cSum <- eventSummary(hm, measures = F)
cSum <- rename(cSum, StartTime = start, EndTime = end)
cSum <- select(cSum, StartTime, EndTime, species, nClicks, Latitude, Longitude, db)
cSum$Project <- 'CCES'
cSum$Deployment <- as.numeric(gsub('.*Drift-([0-9]{1,2}).*', '\\1', basename(cSum$db)))
cSum <- filter(cSum, species %in% c('BB', 'BW', 'BW37V', 'BW43', 'BWC', 'MS', 'ZC', 'BWposs'))
