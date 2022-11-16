source('getClipFromDf.R')
library(tuneR)
# this needs "start" and "end" as POSIXct and "id" to use for labeling
df <- read.csv('ShipDet.csv', stringsAsFactors = FALSE)
# change to whatever format your dates are in
df$start <- as.POSIXct(df$start, format='', tz='UTC')
df$end <- as.POSIXct(df$end, format='', tz='UTC')
# this just gets the clips as WaveMC class objects, point to folder of recordings
clips <- getClipFromDf(df, rec='../Data/3D2PAMr_test_files/Recordings/')

# you can also provide a function that does something that will be
# applied to these clips instead
# "wav" is the WaveMC, "name" is from your "id" column, "time" is
# vector of (start, end), "channel" is channel #, "mode" is not
# relevant but left over from the original "getClipData"
# Any other named arguments can also be provided to "getClipFromDf"
writeClip <- function(wav, name, time, channel, mode, dir='.', ...) {
    if(!dir.exists(dir)) {
        dir.create(dir)
    }
    # you could change this to include the start time of the recording
    # by adding as.character(time[1]) or something like that to the file
    # name here if you want
    file <- file.path(dir, paste0(name, '.wav'))
    writeWave(wav, file, extensible = FALSE)
    # return the file name for keeping track of
    file
}

# this will write clips and put them in folder "dir"
clipFiles <- getClipFromDf(df, '../Data/3D2PAMr_test_files/Recordings/', FUN=writeClip, dir='./TestClips')

# if you want to create spectrograms you can do that in another function (left as exercise for the reader ;) )
