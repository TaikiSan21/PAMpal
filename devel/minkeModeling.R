# Minke whale random forest model dataset creation
# Minke calls have been validated using spec anno module
# in PAMguard, this script reads in all data from validated
# wav files and marks GPL detections that are within annotation boxes

# 2023-04-10: First version TSakai

#### CREATE TRAINING DATA ####
library(PAMpal)
# you can do this with multiple databases and binary folders at once
db <- '../Data/MinkeVal/CCES2018_018_PG_2_02_02_Minke.sqlite3'
bin <- '../Data/MinkeVal/BINARY/'
pps <- PAMpalSettings(db, bin, sr_hz='auto', filterfrom_khz=0, filterto_khz=NULL, winLen_sec=.0025)
# validated files have been added to database, so process with mode='db'
validated <- processPgDetections(pps, mode='db', id='CCES018_Minke')
# this function marks detections with TRUE/FALSE if they are within spec anno boxes
# adds columns "inAnno" (the T/F marker) and "annoId" records the annotation box(es)
# that it matched. annoId is mostly for debugging/validity checking purposes
validated <- markAnnotated(validated)

valGpl <- getGPLData(validated)
# check how many in/out of annotations
table(valGpl$inAnno)

# assign species values based on annotations
valGpl$species <- 'NotMinke'
valGpl$species[valGpl$inAnno] <- 'Minke'
valGpl$species <- factor(valGpl$species)

#### TRAIN MODEL ####
# columns that should not be used for training model
nonModelCols <- c('UTC', 'UID', 'db', 'inAnno', 'annoId', 'detectorName', 'BinaryFile', 'eventLabel', 'eventId')
# build rfPermute model - not sure if this is what we want to do + will need some parameter tuning
library(rfPermute)
# the "setdiff" drops the above columns from the dataset so tehy arent used
rf <- rfPermute(species ~ ., data=valGpl[setdiff(colnames(valGpl), nonModelCols)])
summary(rf)

#### CREATE PREDICTION DATA ####
# process all data with mode='recording' so we can get only unvalidated data
all <- processPgDetections(pps, mode='recording')
valWavs <- sapply(events(validated), function(x) {
    gsub('Added by PAMpal, event ID: ', '', ancillary(x)$eventComment)
})
allWavs <- names(events(all))
notValidated <- all[setdiff(allWavs, valWavs)]
# check nothing weird happened, this should be TRUE
length(events(all)) == (length(events(validated)) + length(events(notValidated)))
predGpl <- getGPLData(notValidated)
predGpl$predProb <- predict(rf, newdata=predGpl[setdiff(colnames(predGpl), nonModelCols)], type='prob')[, 'Minke']
predGpl$predClass <- predict(rf, newData=predGpl[setdiff(colnames(predGpl), nonModelCols)])
