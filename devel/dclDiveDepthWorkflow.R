# Updated 5-20-2024 First revised edition for DCLDE Workshop

source('dclDiveDepthFunctions.R')

#### PAMpal processing ####

# paths to PAMGuard data
db <- '../Data/DCLDE_DD/DiveDepthDetection_Annotated_examples_HB1603_SPWH_Leg2_20160725.sqlite3'
binFolder <- '../Data/DCLDE_DD/20160725/'
# "sr_hz" should be changed to decimated sample rate
pps <- PAMpalSettings(db=db,
                      binaries = binFolder,
                      sr_hz=96e3,
                      filterfrom_khz=2,
                      filterto_khz=NULL,
                      winLen_sec=.0025)
data <- processPgDetections(pps, mode='db')

#### Add required metadata to AcousticStudy ####
# "depthCol" is name of column to use in "Hydrophone_Depth_Data" table, 
# can also be provided as a dataframe
data <- data %>% 
    addGps() %>% 
    addHydrophoneDepth(depthCol='Sensor_0_Depth')
data <- addRecordings(data, folder='../Data/DCLDE_DD/20160725_Soundfiles/')

# save here so we don't have to reprocess in future
saveRDS(data, file='../Data/DCLDE_DD/20160725_Study.rds')

#### Create wave clips ####

# this is folder to store outputs
wavDir <- '../Data/DCLDE_DD/0725WavClips'    
# buffer and channel may be changed, other settings are standard
writeEventClips(data, 
                buffer=c(0, .2),
                channel=5,
                useSample = TRUE,
                mode='detection',
                outDir=wavDir) 

# wave height is used for error estimation, can be provided either as
# single value, timestamped wave height, or timestamped beaufort
data <- addWaveHeight(data, .5)

#### Calculate estimated dive depths ####

# these are species-specific parameters to use for depth estimation process
# freqLow/high - ranges for filter to apply to signal
# minDepth - minimum detection depth for species
# maxRange - maximum detection range for species
# minDepth/maxRange are used (with array depth) to calculate min/max bounds
# for acceptable slant delay times, these can alternatively be directly
# provided as minTime/maxTime
ddParams <- list(
    freqLow = 2e3,
    freqHigh = 16e3,
    minDepth=400,
    maxRange=6500,
    minTime=NULL,
    maxTime=NULL)

# soundSpeed (m/s) and hpDepthError (m) can be adjusted
# nPlot is number of click waveforms to plot, nCol is number of columns for that plot
# outDir will store the summary figure plots for each event and the wav clip plots to review
data <- calcDiveDepth(data,
                      wavDir, 
                      params=ddParams,
                      soundSpeed=1500,
                      hpDepthError = 1,
                      clipLen=.03,
                      nPlot=50,
                      nCol=5,
                      outDir='../Data/DCLDE_DD/0725Outputs192') #92 seconds @0, 120s@50 

#### Review and filter results ####

# filter clicks by a few criteria - maximum estimated depth, minimum autocorrelation magnitude,
# and swim speed. 
# For swim speed, a maximum swim speed "speed" (m/s) is used to eliminate
# impossible clicks if their depth is more than "time" * "speed" apart
data <- filterClicks(data, time=30, speed=50/30, maxDepth=4e3, minCorr = .012) # 1.33s

# launch shiny app to review detections by looking at estimated depth tracks and echograms
data <- runDepthReview(data)

# use the "keepClick" column to filter our data down to only selected, then we can summarise
clicks <- filter(data, keepClick) %>% getClickData

