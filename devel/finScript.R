######################## Processing Section - Can Skip ##########################
library(PAMpal)
db <- '../Data/CHW_Fin/Databases/'
bin <- '../Data/CHW_Fin/Binaries/'
pps <- PAMpalSettings(db, bin, sr_hz=200, filterfrom_khz=0, filterto_khz=0, winLen_sec=1,threshold=-5,
                      functions = list('ClickDetector' = finClick), verbose=FALSE)
data <- processPgDetections(pps, mode='db', id='FinAllValidated')
saveRDS(data, file='../Data/CHW_Fin/FinAllValidated.rds')
sub20csv <- list.files(bin, recursive = TRUE, pattern='20Percent.csv$', full.names=TRUE)
csvDf <- lapply(sub20csv, function(x) read.csv(x, stringsAsFactors = FALSE))
names(csvDf) <- basename(sub20csv)
View(csvDf)

lapply(csvDf, function(x) x$pg_event)

dbNames <- gsub('\\.sqlite3$', '', basename(files(data)$db))
result <- vector('list', length=length(csvDf))
for(i in seq_along(result)) {
    result[[i]] <- list(db=dbNames[i], evNum = csvDf[[i]]$pg_event, evName=paste0(dbNames[i], '.OE', csvDf[[i]]$pg_event))
}
result <- bind_rows(result)
saveRDS(result, file = '../Data/CHW_Fin/20SubDf.rds')
# fuck, this is ONLY Fin+ events/detections. For Fin- we need to read all.
# and then re/un-fuck the event names
goods <- data[unlist(result)]

saveRDS(goods, file = '../Data/CHW_Fin/Fin20Validated.rds')

dbBin <- data.frame(db=list.files('../Data/CHW_Fin/Databases/', full.names=TRUE),
                    bin=list.dirs('../Data/CHW_Fin/Binaries/', recursive = FALSE, full.names = TRUE))

result <- vector('list', length=nrow(dbBin))
for(i in 1:nrow(dbBin)) {
    cat('Processing', i, '\n')
    thisPps <- PAMpalSettings(db=dbBin$db[i], binaries = dbBin$bin[i], sr_hz=200, filterfrom_khz=0, filterto_khz=NULL,
                          winLen_sec=1, threshold=-5,
                          functions = list('ClickDetector' = finClick), verbose=FALSE)
    thisData <- processPgDetections(thisPps, mode='recording', id=paste0(basename(dbBin$bin[i]), '_recording'))
    saveRDS(thisData, file=paste0('../Data/CHW_Fin/', basename(dbBin$bin[i]), 'AllRaw.rds'))
    result[[i]] <- thisData[csvDf[[i]]$pg_event]
}


allRaws <- list.files('../Data/CHW_Fin/', pattern='[0-9]+AllRaw.rds', full.names = TRUE, recursive = FALSE)
result <- vector('list', length=length(allRaws))
csvList <- split(csv20, csv20$db)
for(i in seq_along(allRaws)) {
    thisDrift <- readRDS(allRaws[i])
    has1 <- sapply(events(thisDrift), function(x) {
        any(grepl('1', names(detectors(x))))
    })
    thisDrift <- thisDrift[has1]
    thisDrift <- thisDrift[sapply(events(thisDrift), function(x) nClicks(x) > 0)]
    # weird issue where same wav file was in 3 tiems on drift 21
    if(i == 9) {
        result[[i]] <- thisDrift
        next
    }
    result[[i]] <- thisDrift[csvList[[i]]$evNum]
}
saveRDS(result, '../Data/CHW_Fin/AcList.rds')
all20 <- readRDS('../Data/CHW_Fin/AcList.rds')
all20 <- bindStudies(all20)
saveRDS(all20, '../Data/CHW_Fin/All20.rds')
val20 <- readRDS('../Data/CHW_Fin/FinAllValidated.rds')
csv20 <- readRDS('../Data/CHW_Fin/20SubDf.rds')
val20 <- val20[csv20$evName]
saveRDS(val20, '../Data/CHW_Fin/Val20.rds')

# Read in saved data
val20 <- readRDS('../Data/CHW_Fin/Val20.rds')
val20Wavs <- unname(sapply(events(val20), function(x) gsub('Added by PAMpal, event ID: ', '', ancillary(x)$eventComment)))
all20 <- readRDS('../Data/CHW_Fin/All20.rds')
# SOrting out species IDs - events that had any confirmed fin click are FinWhale
spMark <- rep('NotFin', length(events(all20)))
spMark[which(names(events(all20)) %in% val20Wavs)] <- 'FinWhale'
all20 <- setSpecies(all20, method='manual', value=spMark)
# Calculate ICI for these and build banter model
all20 <- calculateICI(all20)
# Adding detector level ici (time since last call, both within detector and overall)
for(e in seq_along(events(all20))) {
    thisEv <- all20[[e]]
    if(nClicks(thisEv) == 0) {
        next
    }
    ici <- getICI(thisEv, 'data')
    for(d in seq_along(detectors(thisEv))) {
        if(!names(detectors(thisEv))[d] %in% ici$detectorName) {
            next
        }
        detIci <- filter(ici, detectorName == names(detectors(thisEv))[d])
        detIci$detectorIci <- detIci$ici
        allIci <- filter(ici, detectorName == 'All')
        allIci$allIci <- allIci$ici
        detectors(thisEv)[[d]] <- left_join(detectors(thisEv)[[d]], detIci[c('UID', 'detectorIci')], by='UID')
        detectors(thisEv)[[d]] <- left_join(detectors(thisEv)[[d]], allIci[c('UID', 'allIci')], by='UID')
    }
    all20[[e]] <- thisEv
}

saveRDS(all20, '../Data/CHW_Fin/All20Final.rds')

###################### Processing Section Done ########################

all20 <- readRDS('../Data/CHW_Fin/All20Final.rds')
bntData <- export_banter(all20)
# Drop WMD
bntData$detectors$Whistle_and_Moan_Detector <- NULL
library(banter)
bnt <- initBanterModel(bntData$events)
bnt <- addBanterDetector(bnt, bntData$detectors, ntree=500, sampsize = 50)
bnt <- runBanterModel(bnt, ntree=500, sampsize = 50)
summary(bnt)
rfPermute::plotVotes(bnt@model)
# Getting True Pos / etc. AcousticStudies

val20 <- readRDS('../Data/CHW_Fin/Val20.rds')
# Validated & Detector1 is True Positive
truePos <- filter(val20, grepl('1$', detector))
# Validated & Detector0 is False Negative
falseNeg <- filter(val20, grepl('0$', detector))
# Detector 1 & Not In Validated is False Positive
falsePos <- filter(all20, grepl('1$', detector))
falsePos <- filter(falsePos, !UID %in% getClickData(truePos)$UID)
# Detector 0 & Not In Validated is True Negative
trueNeg <- filter(all20, grepl('0$', detector))
trueNeg <- filter(trueNeg, !UID %in% getClickData(falseNeg)$UID)
# these are close to total but slightly off? 6 less here than in total
# there are overlapping UIDs or some are marked as both 0/1 i think
nClicks(truePos)
nClicks(falseNeg)
nClicks(trueNeg)
nClicks(falsePos)

# Plots for fins - need to update to latest GitHub PAMpal for one of these options to work
png('TruePositiveAllDriftsConcat.png', width=6, height=4, units='in', res=300)
calculateAverageSpectra(truePos, 1:10e3, wl=128, noise=TRUE, showBreaks = FALSE, plot=c(T,F),
                        title='True Positive All Drifts')
dev.off()
png('TruePositiveAllDriftsAvgSpec.png', width=6, height=4, units='in', res=300)
calculateAverageSpectra(truePos, 1:10e3, wl=128, noise=TRUE, showBreaks = FALSE, plot=c(F,T),
                        title='True Positive All Drifts')
dev.off()
png('FalseNegAllDriftsConcat.png', width=6, height=4, units='in', res=300)
calculateAverageSpectra(falseNeg, 1:10e3, wl=128, noise=TRUE, showBreaks = FALSE, plot=c(T,F),
                        title='False Negative All Drifts')
dev.off()
png('FalseNegAllDriftsAvgSpec.png', width=6, height=4, units='in', res=300)
calculateAverageSpectra(falseNeg, 1:10e3, wl=128, noise=TRUE, showBreaks = FALSE, plot=c(F,T),
                        title='False Negative All Drifts')
dev.off()
# plots for not fins
png('FalsePosAllDriftsConcat.png', width=6, height=4, units='in', res=300)
calculateAverageSpectra(falsePos, 1:10e3, wl=128, noise=TRUE, showBreaks = FALSE, plot=c(T,F),
                        title='False Positive All Drifts')
dev.off()
png('FalsePosAllDriftsAvgSpec.png', width=6, height=4, units='in', res=300)
calculateAverageSpectra(falsePos, 1:10e3, wl=128, noise=TRUE, showBreaks = FALSE, plot=c(F,T),
                        title='False Positive All Drifts')
dev.off()
png('TrueNegAllDriftsConcat.png', width=6, height=4, units='in', res=300)
calculateAverageSpectra(trueNeg, 1:10e3, wl=128, noise=TRUE, showBreaks = FALSE, plot=c(T,F),
                        title='True Negative All Drifts')
dev.off()
png('TrueNegAllDriftsAvgSpec.png', width=6, height=4, units='in', res=300)
calculateAverageSpectra(trueNeg, 1:10e3, wl=128, noise=TRUE, showBreaks = FALSE, plot=c(F,T),
                        title='True Negative All Drifts')
dev.off()
