# Updated 2024-11-12: Adjusting for combo whistle+click export example
# Updated 2024-10-07: ICI for bins now implemented + filtering example
source('soundscopeFunctions.R')

# Process with PAMpal ####

db <- '../Data/DCLDE_DD/DiveDepthDetection_Anntotated_examples_HB1603_SPWH_Leg1_20160707.sqlite3'
bin <- '../Data/DCLDE_DD/20160707/'
# these two can be however you normally process the data
pps <- PAMpalSettings(db=db, binaries = bin, sr_hz=96e3, filterfrom_khz=2, filterto_khz=NULL,
                      winLen_sec=.0025)
data <- processPgDetections(pps, mode='recording')
# this step points to recording folders, just like with the dive depth workflow
rec <- '../Data/DCLDE_DD/20160707_Soundfiles/'
data <- addRecordings(data, folder=rec)
# Example of save/load so you don't have to re-process
# data <- saveRDS('../Data/DCLDE_DD/20160707_Study_AllRecording.rds')
# data <- readRDS('../Data/DCLDE_DD/20160707_Study_AllRecording.rds')

# Export Clicks ####
# bin - sets how long of a time period to use for binning click detections, and
#   this is how long the spectrogram will be in soundscope. I think max is 
#   20 seconds, but I'm not sure
# binIci - if TRUE and "bin" is not NULL calculates the ICI of each time bin
# extraCols - these aren't yet used by Soundscope, so don't worry about them. But
#   eventually these will show up as parameters you can view alongside the 
#   spectrogram. You can name any of the PAMpal calculated values here.
scopeClicks <- export_soundscope(data, 
                               detector='click', 
                               bin='10second',
                               binIci=TRUE,
                               extraCols=c('peak', 'dBPP'))

createSoundscopeNc(data=scopeClicks, file='testSoundscopeClicks.nc')

# Export Whistles ####
# duration and frequency_bandwidth should be default exported measures 
# so no "extraCols" needed
scopeWhistles <- export_soundscope(data,
                                   detector='whistle')

createSoundscopeNc(data=scopeWhistles, file='testSoundscopeWhistles.nc')

#### Example to filter out certain binIci values ####
# Probably best to plot the ICIs and look at them before trying to filter
library(ggplot2)
ggplot(scopeClicks$data, aes(x=n, y=binIci)) +
    geom_point()
scopeClicks$data <- filter(scopeClicks$data, binIci > 0.1)
createSoundscopeNc(data=scopeClicks, file='testSoundscopeClicksIciFilter.nc')


#### Fixing recording paths of an existing nc file ####
ncFile <- '../Data/Soundscope/NEFSC_SBNMS_202208_SB01_10secBin_Soundscope.nc'
# point to folder of recordings with correct path
recFolder <- 'Y:/BOTTOM_MOUNTED/NEFSC_SBNMS/NEFSC_SBNMS_202208_SB01/67403784_48kHz_UTC'
# this should tell you how many it fixed
renameSoundscopeRecordings(ncFile, recFolder)
