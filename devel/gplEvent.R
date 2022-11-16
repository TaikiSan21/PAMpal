db <- '../Data/GPL_Event/GPL_BlueAB_PG2_02_03_CCES2018_Drift21_ClickDetectorAdded.sqlite3'
bin <- '../Data/GPL_Event/Drift_21/'
library(PAMpal)
pps <- PAMpalSettings(db, bin, sr_hz='auto', filterfrom_khz=0, filterto_khz=NULL, winLen_sec=.0025)
# data <- processPgDetections(pps, mode='recording')
# saveRDS(data, '../Data/GPL_Event/GPLStudy.rds')
data <- readRDS('../Data/GPL_Event/GPLStudy.rds')
nDets <- sapply(events(data), nGPL)
hasDets <- data[nDets > 0]

set.seed(1296)
sub20 <- sample(names(events(hasDets)), .2*length(events(hasDets)), replace=FALSE)
library(PAMmisc)

for(e in seq_along(events(hasDets))) {
# for(e in 1) {
    thisEv <- hasDets[[e]]

    addPgEvent(db=db, binary = files(thisEv)$binaries, UIDs = getGPLData(thisEv)$UID,
               eventType = ifelse(id(thisEv) %in% sub20, 'Validate', 'Ignore'),
               comment = paste0('Created by PAMpal, event ID: ', id(thisEv)))
}

df <- data.frame(pamguardId = 1+1:length(events(hasDets)),
                 wavFile = names(events(hasDets)),
                 eventType = ifelse(names(events(hasDets)) %in% sub20, 'Validate', 'Ignore'))
df$done <- ifelse(df$eventType == 'Validate', '', 'x')

write.csv(df, file='../Data/GPL_Event/Drift21_20pctValidate.csv', row.names = FALSE)
