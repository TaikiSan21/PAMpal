source('devel/matchTemplateFunctions.R')
# change if the set of templates changes
templateNames <- c("ZC","BW43","BW39V","MS","BB","BW70")
# keeping the match/reject values just in case they are useful, but
# can remove these in future if they aren't necessary
extraCols <- c(paste0(templateNames, '_match'),
               paste0(templateNames, '_reject'))

baseDir <- '../Data/MTC_CV/Binaries/ADRIFT_024/'
binFolder <- baseDir
# this database should be a COPY of the original because we will add events to it later
db <- '../Data/MTC_CV/Databases/ADRIFT_024_MTC.sqlite3'
# the binary processing takes a really long time, this automatically saves to an RDS file
# so that you don't have to reprocess in future
saveFile <- file.path(baseDir, 'Drift24Template.rds')

allData <- loadTemplateFolder(binFolder, names=templateNames, extraCols=extraCols, file=saveFile)
# allData <- readRDS(saveFile)
# these are in order of "templateNames" above. Can look at data and see if any of these need to
# be raised/lowered
# threshVals <- c(.06, .15, .15, .15, .15, .15)
threshVals <- c(.55, .45, .6, .5, .25, .6)
# change "db" to "db=NULL" if we don't have manual events that we can match UIDs to

allData <- addTemplateLabels(allData, db, templateNames, threshVals)
# minDets is minimum detections to count as an event, maxSep is max time between detections
# before an event is ended. maxLength is maximum length of an event
allData <- markGoodEvents(allData, minDets=3, maxSep=120, maxLength=120)

# summary of how many of the detections in manually annotated events were tagged by template
# Adding the "db" argument will include the time-overlap based summary
manualSummary <- summariseManualEvents(allData, db=db)
# summary of how many detections tagged by template were present in manually annotated events
# Adding the "db" argument will include the time-overlap based summary
templateSummary <- summariseTemplateEvents(allData, db=db)

# Labeling "allData" by the template event results
# Don't want to use all FP template events (these are labeled "none"),
# subset randomly some fixed amount instead
whichNotBW <- which(templateSummary$overlapLabel == 'none')
# good practice to set seed before any randomness for reproducibility
set.seed(112188)
# change this to however many NotBW events you want
nSubset <- 200
whichSubset <- sample(whichNotBW, size=nSubset, replace=FALSE)
# we'll only label the selected subset as "NotBW", leave the rest as "none"
templateSummary$overlapLabel[whichSubset] <- 'NotBW'
# attach these template labels to the data
allData <- left_join(allData, templateSummary[c('templateEvent', 'overlapLabel')], by='templateEvent')
# remove all "none" events so we only add the labeled and selected subset to database as events
addTemplateEvents(db, binFolder, filter(allData, overlapLabel != 'none'), labelCol='overlapLabel')

# adds events meeting nDets/nSeconds criteria to the database
# make sure db is a COPY of the original for safety
# commented out in favor of above approach
addTemplateEvents(db, binFolder, allData)

### OPTIONAL process again with PAMpal to do stuff ###
library(PAMpal)
pps <- PAMpalSettings(db, binFolder, sr_hz=288e3, filterfrom_khz=10, filterto_khz=NULL, winLen_sec=.0025)
data <- processPgDetections(pps, mode='db', id='MatchTempAD24')
data <- setSpecies(data, method = 'pamguard')
# new "FP" and "TP" events in addition to originals
table(species(data))
saveRDS(data, 'D:/CCES_2018/Finalized BW Analyses/Drift-19 (completed by JST)/12 dB threshold/matchTemplateStudy19.rds')


#### Template event filtering tests

summaryPlotData <- function(x, templateNames, threshVals) {
    matchCols <- paste0(templateNames, '_match')
    x <- split(x, x$templateEvent)
    x <- bind_rows(lapply(x, function(ev) {
        if(ev$templateEvent[1] == 'none') {
            return(NULL)
        }
        result <- vector('list', length=length(matchCols))
        names(result) <- matchCols
        for(i in seq_along(result)) {
            result[[i]] <- mean(ev[[names(result)[i]]], na.rm=TRUE)
        }
        result$nDets <- nrow(ev)
        result$templateEvent <- ev$templateEvent[1]
        result$templateGood <- ev$templateGood[1]
        evLabel <- unique(ev$eventLabel)
        if(length(evLabel) > 1) {
            evLabel <- evLabel[evLabel != 'none']
        }
        result$eventLabel <- evLabel
        result
    }))
    x <- tidyr::pivot_longer(x, cols=all_of(matchCols), names_to='Template', values_to='MeanScore')
    x
}
allData$templateGood <- allData$templateEvent %in% templateSummary$templateEvent[templateSummary$pctOverlap > 0]
spd <- summaryPlotData(allData, templateNames, threshVals)
ggplot(spd, aes(x=nDets)) +
    geom_point(aes(y=MeanScore, col=eventLabel), size=2) +
    facet_wrap(templateGood~Template, ncol=6) +
    ggtitle('Detections vs Mean Sore')
spd %>%
    select(-Template, -MeanScore) %>%
    distinct() %>%
    ggplot(aes(x=nDets)) +
    geom_histogram(aes(fill=eventLabel)) +
    facet_wrap(~templateGood, scales='free') +
    xlim(0, 20) +
    ggtitle('Number of Detections')

