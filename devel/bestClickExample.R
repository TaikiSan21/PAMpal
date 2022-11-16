# I misspoke about being able to re-use the same parameters, that actually didn't work
# yet in PAMpal. So I changed it, this now works if you update to GitHub version 0.15.1.
# The new feature is being able to specify new functions when you first call `PAMpalSettings`
# With the "functions" argument

# from_khz sets the lower bound of freq range to search
# to_khz sets the upper bound of freq range to search
# Then it just checks which channel has a higher max dB in that range
# If only one channel, always marked as best
pps <- PAMpalSettings(db='../Data/AMTask/3D2PAMr_test_files/HB1603_MF_Master_Leg3_Skala.sqlite3',
                      binaries = '../Data/AMTask/3D2PAMr_test_files/BinaryFiles_20160815_WithUID/',
                      sr_hz = 'auto',
                      filterfrom_khz=10,
                      winLen_sec=.0025,
                      filterto_khz=NULL,
                      functions=list('ClickDetector' = markBestClick),
                      from_khz=30,
                      to_khz=50)

# Of if you dont want to update, remove last 3 lines above (from functions and below) and add with addFunction
# If doing it this way you unfortunately have to specify these parameters again
pps <- addFunction(pps, markBestClick, 'ClickDetector', from_khz=30, to_khz=50, sr_hz='auto', winLen_sec=.0025)

# Then process data like normal
data <- processPgDetections(pps)

# This function just adds a new column "bestClick" that is TRUE if best click FALSE if not.
# If only one channel was present, it is marked TRUE.

bestClicksOnly <- filter(data, bestClick == TRUE)

# It doesn't really matter since they will all be TRUE (so no information), but probably
# want to remove this variable from any BANTER data:
banter <- export_banter(bestClicksOnly, dropVars = c('bestClick'))