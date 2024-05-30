source('dclDiveDepthFunctions.R')
# paths to pamguard data
db <- '../Data/diveDepth/HB2102_Leg1_06262021_SPWH_Database.sqlite3'
bin <- '../Data/diveDepth/Binaries/20210626/'
# PAMpal processing - "sr_hz" will need to be changed to decimated sample rate
pps <- PAMpalSettings(db, bin, sr_hz=96e3, filterfrom_khz=2, winLen_sec=.0025, filterto_khz=NULL)
data <- processPgDetections(pps,mode='db')
data <- addGps(data)
# "depthCol" is column name in "Hydrophone_Depth_Data" table
# can also be provided as dataframe
data <- addHydrophoneDepth(data, depthCol='Sensor_0_Depth')
data <- addRecordings(data,
                      folder='Z:/Stellwagen_Old/ACOUSTIC_DATA/TOWED_ARRAY/NEFSC_NE-OFFSHORE_202106_HB2102/Leg1/20210626')
# save here so that we don't need to reprocess
saveRDS(data, file='../Data/diveDepth/AcSt20210626.rds')
data <- readRDS('../Data/diveDepth/AcSt20210626.rds')
wavDir <- '../Data/diveDepth/20210626Wavs/'
# wavs <- writeEventClips(data,
#                         buffer = c(0, 0.2),
#                         outDir = wavDir,
#                         mode = 'detection',
#                         channel = 3,
#                         useSample = TRUE)


dd <- calcDiveDepth(data, wavDir, 
                    params=list(
                        freqLow = 2e3,
                        freqHigh = 16e3,
                        minDepth=400,
                        maxRange=6500),
                    soundSpeed=1500,
                    nCol=5,
                    nPlot=0,
                    hpDepthError = 1,
                    clipLen=.03,
                    outDir='../Data/diveDepth/20210626Test')

saveRDS(dd, file='../Data/diveDepth/AcSt20210626_wDD.rds')
dd <- readRDS('../Data/diveDepth/AcSt20210626_wDD.rds')
dd <- filterClicks(dd, time=30, speed=50/30, maxDepth=4000, minCorr=.012)
dd <- runDepthReview(dd)


