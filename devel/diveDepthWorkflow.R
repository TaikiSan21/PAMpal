####### CHANGE LOG #######
# 05/19/2021 Dive Depth Automation - Needs PAMpal v 0.12.5
# 05/24/2021 Updated to remove events with NA localizations @ line 31-40
# 6/21/2022 Wow its been a year. Updating to be more flexible with wav height adding step
##########################

# change this to location of R file
source('./devel/diveDepthFunctions.R')
# DB and Binaries for processing
db <- '../Data/AMTask/3D2PAMr_test_files/HB1603_MF_Master_Leg3_Skala.sqlite3'
bin <- '../Data/AMTask/3D2PAMr_test_files/BinaryFiles_20160815_WithUID/'

# Run PAMpal (any settings here can be changed to your desired values)
pps <- PAMpalSettings(db=db, binaries=bin, sr_hz='auto', filterfrom_khz=5, filterto_khz=NULL, winLen_sec=.0025)
# id here will be used for metadata csv name
data <- processPgDetections(pps, mode='db', id='AMTest')

# Add GPS, Depth, and recording folder info
data <- addGps(data) %>%
  # addHydrophoneDepth() %>%
  addRecordings(folder = '../Data/3D2PAMr_test_files/Recordings/', log=FALSE)
# this will by default look for "Hydrophone_Depth_Data" table in database, and use column "depthCol" for the data
# Can also provide as a df to "depth", needs column "UTC". Matches to nearest, may need to change "thresh" value
data <- addHydrophoneDepth(data, depthCol = 'Sensor_0_Depth')
# FIltering down to only BEAK, POBK, and PRBK events
# Since these are currently in eventType, easiest way is to temporarily assign
# them as the species ID and use PAMpal's filtering function, then we
# reset the species IDs back to NA
goodData <- data %>%
  setSpecies(method='pamguard') %>%
  filter(species %in% c('BEAK', 'POBK', 'PRBK')) %>%
  setSpecies(method='manual', value=NA)

# Checking for events that dont have NAs for their TargetMotion locs
hasLocalization <- sapply(events(goodData), function(x) {
  tarMo <- localizations(x)$PGTargetMotion
  if(is.null(tarMo)) {
    return(FALSE)
  }
  !is.na(tarMo$TMLatitude1) &
    !is.na(tarMo$TMPerpendicularDistance1)
})
goodData <- goodData[hasLocalization]

# Now for the real species IDs, this will attempt to parse the comments and write
# a CSV file. It will tell you how many it could not find an easy match for,
# and those will be left as NAs
commSpec <- commentToSpecies(goodData, validate = 'SpeciesValidation.csv')
# Now open up the "SpeciesValidation.csv", check all NAs to see if they should
# be something else. Make a note of any cases that get missed, I can add them to
# the special cases the function checks.

# When done with that, read it back in
commSpec <- read.csv('SpeciesValidation.csv', stringsAsFactors = FALSE)
# We assign these species manually, then remove any that were left as NA.
#
goodData <- setSpecies(goodData, method='manual', value=commSpec) %>%
  filter(!is.na(species))

# Folder where wav clips should be stored, if it doesnt exist yet funciton will create it
clipDir <- '../Data/3D2PAMr_test_files/WavClips'
wavs <- writeEventClips(goodData,
                        buffer = c(0, 0.2),
                        outDir = clipDir,
                        mode = 'detection',
                        channel = 5,
                        useSample = TRUE)

# effort can be .xlsx or .csv of just that sheet
effortTable <- '../Data/3D2PAMr_test_files/HB1603beaufort_effort.xlsx'
# This will create the metadata CSV table and the PNG of waveform clips
# It will also store both of these outputs - the csv as a dataframe, and the wav clips
# as a matrix with 1 column for each wavform. The wave stored is the exact shortened
# and filtered one shown in the PNG - maybe these could be passed directly to other code?

# PARAMETERS:
# "outDir" is where the CSV and PNG will be created, if NULL will default to the wavFolder
# "file" is the name of the CSV file, if NULL will default to "DiveDepthMetadata_STUDYID.csv"
#    where STUDYID is the id set in line 8 aobve (AMTest in this example)
# "wavFolder" folder where wav clips can be found
# "effort" effort file location
# "clipLength" length of clip for PNG. Can be seconds or samples, it will sort it
#    out based on the size of the input (if < 10 -> assume seconds, else samples)
# "soundSpeed" in m/s
# "depthSensAcc" accuracy term in table

height <- readBftNEFSC(effortTable)
depthOutputs <- export_diveDepthNEFSC(goodData,
                              outDir = NULL,
                              file = NULL,
                              wavFolder = clipDir,
                              height = height,
                              clipLength = .03,
                              soundSpeed = 1500,
                              depthSensAcc = 1)

