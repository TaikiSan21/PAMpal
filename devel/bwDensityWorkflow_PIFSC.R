# Beaked whale density workflow for PIFSC

library(PAMpal)
source('../Data/RoboJ/Taikis FIles/RoboJayFunctions.R')
# Place to save outputs
outPath <- '../Data/RoboJ/Test'
# Can test on just one station, if that works then we can run all stations at once
db <- '../Data/RoboJ/Test/Station-7_Card-A_MASTER-BW - Copy.sqlite3'
bin <- '../Data/RoboJ/Test/Binaries_Station-7_Card-A_WithUID/'

pps <- PAMpalSettings(db, bin, sr_hz='auto', winLen_sec=.0025, filterfrom_khz=10, filterto_khz=NULL)
data <- processPgDetections(pps, mode='db')
data <- setSpecies(data, method='pamguard')
# Re-label species codes if necessary for consistency
speciesRename <- data.frame(old = c('BW34-50', 'BW46', 'BW50-75', '?BW', 'Zc', 'Pm', 'SW', '?Pm', 'Oo'),
                            new = c('BW37V', 'MS', 'MS', 'BWunid', 'ZC', 'PM', 'PM', '?PM', 'OO'))
data <- setSpecies(data, method = 'reassign', value=speciesRename)
# Make sure this doesnt have any codes that didnt get converted above
unique(species(data))

data <- addGps(data)
# Folder containing wavs. Easy if just one DB:
recordingFolder <-  '../Data/RoboJ/Test/Recordings'
# If multiple DBs, first check order of DBs:
files(data)$db
# Then folder needs to be a vector in that order specifying which folder for each DB:
recordingFolder <- c(
    '../Data/RoboJ/Test/RecordingsDB1',
    '../Data/RoboJ/Test/RecordingsDB2'
)
# Alternatively set folder=NULL and it will do a pop-up asking you to choose folder for each DB if thats easier
data <- addRecordings(data, folder=recordingFolder)
# saveRDS(data, file=file.path(outPath, 'RoboJayStudy.rds'))

# stationPattern is a REGEXP to turn event names into a single "Station" name that is
# easier to write down.
# Default gets the name of the database from PAMpal's eventIds that are created when processing mode='db'
defaultPat <- '(.*)\\.OE.*|\\.DGL.*'
# If your database names contain a single number that labels each drift uniquely, that is probably
# an easier way to label them. Here's an example of the pattern used to extract that number
# from DBs from the CCES survey
CCESPat <- '.*Drift-([0-9]{1,2})_.*[\\.OE|\\.DGL][0-9]+$'
eventClicks <- formatClickDf(data, stationPattern = defaultPat)

# Preps for the soundspeed profile download - marks detections into different groups
# first by the "splitBy" column, then new groups if separated by longer than "time"
# If download times seem to be long you can increase the "time" value to maybe 3 or 7 days
eventClicks <- splitSSPGroups(eventClicks, splitBy='station', time=3600*24*1)

# This creates a SSP for each row, has to download from HYCOM server so may be slow
# Will save the downloaded data to "file" so that it can be read back in and used in the future
sspList <- createSSPList(eventClicks, file=file.path(outPath, paste0(as.character(Sys.Date()), '_sspList.rds')))

# Can read this from disk if running again in future on same data (change date for future)
# sspList <- readRDS(file.path(outPath, '2022-04-20_sspList.rds'))

# change this to value appropriate for your survey
hpDepth <- 110
# these are values for Ziphius from Jay, will need to change for diff species
animalDepth <- 1217
sdAnimalDepth <- 354
# Uses the SSP to correct angles using raytrace, then applies raytrace angle correction
# to existing angles. May be quite slow since raytracing can take some time
eventClicks <- doAngleCorrection(eventClicks, sspList, hpDepth = hpDepth, animalDepth = animalDepth)

# "recorderInfo" needs 5 columns:
# "station" must exactly match the "station" in eventClicks.
# "dutyCycle" is 'minutes on/minutes off', so if recording 2 out of every 10 minutes '2/8'
# "recorder" is just the name of the recorder instrument type
# "deployTime" / "retrTime" are deployment and retrieval times (in UTC). These should be
recorderInfo <- read.csv('../Data/RoboJ/Test/recorderInfo.csv', stringsAsFactors = FALSE)
# change the format here to match whatever your date format is
recorderInfo$deployTime <- as.POSIXct(recorderInfo$deployTime, format='%Y-%m-%d  %H:%M:%S', tz='UTC')
recorderInfo$retrTime <- as.POSIXct(recorderInfo$retrTime, format='%Y-%m-%d  %H:%M:%S', tz='UTC')

eventSummary <- formatEventSummary(eventClicks, recorderInfo = recorderInfo)

recordingDf <- files(data)$recordings

# This output shouldnt take any time to run
snapshots <- tallySnapshots(recordingDf, eventSummary, recorderInfo)

# Do any desired filtering after this point, examples below, species is most important
eventSummary <- filter(eventSummary, species == 'ZC')
# eventSummary <- eventSummary[eventSummary$Station <= 22, ]

# Also takes no time to run
detHist <- createDetectionHistory(eventSummary)

# This will take quite a long time. To test and make sure things just run you can set doJackknife = FALSE and nSim = 1e6,
# but for full run (> 1hr ) we want doJackknife=TRUE and nSim=1e7
detFunction <- newEstDetFunction(eventSummary,
                                 subsetCol = 'recorder',
                                 meanDepth = animalDepth,
                                 sdDepth = sdAnimalDepth,
                                 hpDepth = hpDepth,
                                 truncDist = 4000, # trunc distance for detection function (meters)
                                 sdAngleErr = 0.59, # not sure how Jay settled on this
                                 depthDistr = 'log-normal',
                                 doJackknife = FALSE, # run first with FALSE to make sure everything works, then after set TRUE
                                 jk_nSamples = 20,
                                 nSim=1e7, # can set this lower to 1e5 for a test run
                                 progress=TRUE)
# save all outputs for next steps
saveRDS(eventSummary, file.path(outPath, 'eventSummary.rds'))
saveRDS(detFunction, file.path(outPath, 'detFunction.rds'))
saveRDS(snapshots, file.path(outPath, 'snapshots.rds'))
saveRDS(detHist, file.path(outPath, 'detHist.rds'))