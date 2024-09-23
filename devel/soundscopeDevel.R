# soundscope testo

ncFile <- '../Data/Soundscope/haddock_detections_NEFSC_SBNMS_200601_CH3_newpath.nc'
library(ncdf4)
nc <- nc_open(ncFile)
View(nc)
nc$dim$date$len
# Notes ####
# only units on any variables are for dates and theyre all minutes/seconds/microseconds/milliseconds
# since YYYY-MM-DD HH:MM:SS.UUUUUU
# the "since" dates all seem arbitrary??
# everything has same size - i think were just awkwardly recreating a dataframe
# where the only dimension column is UTC time ("date")

# need to connect to recording files

db <- '../Data/DCLDE_DD/DiveDepthDetection_Anntotated_examples_HB1603_SPWH_Leg1_20160707.sqlite3'
bin <- '../Data/DCLDE_DD/20160707/'
rec <- '../Data/DCLDE_DD/20160707_Soundfiles/'

library(PAMpal)

pps <- PAMpalSettings(db=db, binaries = bin, sr_hz=96e3, filterfrom_khz=2, default=TRUE)
data <- processPgDetections(pps, mode='recording')
dataDb <- processPgDetections(pps, mode='db')
data <- addRecordings(data, folder=rec)
saveRDS(data, '../Data/DCLDE_DD/20160707_Study_AllRecording.rds')
data <- readRDS('../Data/DCLDE_DD/20160707_Study_AllRecording.rds')
data %>% 
    getClickData() %>% 
    ggplot() +
    geom_density(aes(x=peak))

# currently it is using "duration" as the length of the stored
# wav clip - is that useful?
# currently storing Fs as decimated, this conflicts with wav file
# it will actually be reading from. Seems important.
# handling UTC offset

# oooh we will need to pick a single channel for measurements. 
# db only stores channel-agnostic data

# what do we want to do about different detection types (whistle click etc)
# potentially same millisecond time which maybe borks tho this should be rare

# "measurements_name" is what is allowed to browse? or?


hm <- export_soundscope(data, detector='click', channel=5, bin='10second', extraCols=c('peak', 'dBPP'))
# hm$data$audio_channel <- 2
hm$data$frequency_min <- 100.
hm$data$frequency_max <- 16000.
hm$data$frequency_bandwidth <- hm$data$frequency_max-hm$data$frequency_min
createSoundscopeNc(data=hm$data, attributes = hm$attributes, file='testBin.nc')

click <- getClickData(data)
click <- joinRecordingData(click, recData)
click <- click[!is.na(click$file), ]
click <- click[click$Channel ==5, ] #matches now
?plotGram
ix <- 147
plotGram(data, evNum=click$eventId[ix], 
         start=click$UTC[ix]-5, wl=2048,
         end=click$UTC[ix]+5, sr = 96e3, freqRange=c(.1, 16), detCol='red', 
         channel = as.numeric(click$Channel[ix]))

# trying to compare to known PM events 
click <- getClickData(data)
clickDb <- getClickData(dataDb)
evRanges <- lapply(split(clickDb, clickDb$eventId), function(x) range(x$UTC))
click$inMarked <- FALSE
for(i in seq_along(evRanges)) {
    isIn <- click$UTC <= evRanges[[i]][2] & click$UTC >= evRanges[[i]][1]
    click$inMarked[isIn] <- TRUE
}
click <- filter(click, inMarked)
click$inMarked <- click$UID %in% clickDb$UID
# problem is that not all Pms are actually labeled, so data is bad
ggplot(click, aes(x=UTC, color=inMarked)) +
    geom_point(aes(y=angle))
