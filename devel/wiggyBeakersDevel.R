# try wiggies on old AM data
data <- readRDS('../Data/3D2PAMr_test_files/AMStudy.rds')
data <- updateFiles(data, db='../Data/3D2PAMr_test_files/',
                    bin = '../Data/3D2PAMr_test_files/BinaryFiles_20160815_WithUID/')
calculateAverageSpectra(data, evNum = 1:4, plot=c(T,F))

source('../../beaker-wigner-class/R/CV4EProcessingFunctions.R')

amWig <- processAllWig(data, wl=128, sr=192e3, dir='../Data/3D2PAMr_test_files/Wiggies', filter=10)
data <- calculateICI(data, time='peakTime')
ici <- getICI(data)
ici <- bind_rows(lapply(ici, function(x) x['All_ici']), .id='eventId')
ici <- rename(ici, ici = All_ici)
amWig <- left_join(amWig, ici)
amWig$file <- basename(amWig$file)
amWig$species <- 'ZC'
write.csv(amWig, file = '../Data/3D2PAMr_test_files/WigTest.csv', row.names = FALSE)

bigWig <- which(allWigNp$wigMax < .0000001)
bigBin <- getBinaryData(data, UID=amWig$UID[1])
clip <- downSamp(bigBin[[1]]$wave[, 1], srFrom = bigBin[[1]]$sr, 192e3)
clip <- PAMpal:::clipAroundPeak(clip, 128)
plot(clip, type='l')
wiggy <- wignerTransform(clip, n=128, sr=bigBin[[1]]$sr, plot=T)
plot(wiggy$tfr[1,], type='l')

# subset of new data
# SHIT IT SHOULD BE 288e3 FOR FIRST RUN THEN MAKE IT 192e3 LATER
library(PAMpal)
pps <- PAMpalSettings(db=c('../Data/MTC_CV/Databases/PG_2_02_09_ADRIFT_050 - Copy.sqlite3',
                           '../Data/MTC_CV/Databases/PG_2_02_09_CCES_012 - Copy.sqlite3'),
                      binaries = c('../Data/MTC_CV/Binaries/ADRIFT_050/',
                                   '../Data/MTC_CV/Binaries/PG_2_02_09_CCES_012/'),
                      sr_hz=192e3, filterfrom_khz=10, filterto_khz=NULL, winLen_sec=.0025)
dbs <- list.files('../Data/MTC_CV/Databases/', full.names=TRUE)
bins <- list.dirs('../Data/MTC_CV/Binaries/', full.names=TRUE, recursive = FALSE)
pps <- PAMpalSettings(db=dbs, binaries = bins,
                      sr_hz=288e3, filterfrom_khz=10, filterto_khz=NULL, winLen_sec=.0025)
data <- processPgDetections(pps, mode='db')
data <- calculateICI(data, time='peakTime')
saveRDS(data, file='../Data/MTC_CV/ADRIFTStudy_Train288.rds')
data <- readRDS('../Data/MTC_CV/ADRIFTStudy_Train288.rds')
# data <- setSpecies(data, 'pamguard')
# data <- setSpecies(data, method='reassign', value=data.frame(old='?BW', new='PossBW'))
# spDf <- distinct(getClickData(data)[c('eventId', 'species')])
data <- filter(data, Channel == '1')
adWig <- processAllWig(data, wl=128, sr=192e3, dir='../Data/MTC_CV/Wiggies288', filter=10, channel=1)
# oldCsv <- readr::read_csv('../Data/MTC_CV/ADRIFT_Wigner.csv')
ici <- getICI(data)
ici <- bind_rows(lapply(ici, function(x) x['All_ici']), .id='eventId')
ici <- rename(ici, ici=All_ici)
adWig <- left_join(adWig, ici)
adWig$file <- paste0('FP/', basename(adWig$file))
adWig$species <- NULL
write.csv(adWig, file='../Data/MTC_CV/ADRIFT024_Wigner288.csv', row.names=FALSE)
oldCsv$species <- NULL
oldCsv$UID <- as.character(oldCsv$UID)
oldCsv$Channel <- as.character(oldCsv$Channel)
oldCsv$wigUID <- as.character(oldCsv$wigUID)
combined <- bind_rows(oldCsv, adWig[colnames(oldCsv)])
combined$species <- combined$eventLabel
combined$species[!combined$species %in% c('ZC', 'MS', 'BW43', 'BB', 'BW37V')] <- 'BW37V'

write.csv(combined, file='../Data/MTC_CV/ADRIFT_Wigner.csv', row.names = FALSE)
# adWig <- left_join(adWig, spDf)
# adWig$species[adWig$species != 'ZC'] <- 'BW43'
write.csv(adWig, file='../Data/MTC_CV/ADRIFT_Wigner.csv', row.names = FALSE)

preds <- readr::read_csv('../../beaker-wigner-class/preds/adrift_test/adrift_test/ADRIFT_Wigner_pred.csv')

spDf <- distinct(getClickData(data)[c('eventId', 'species')])
preds$species <- NULL
preds <- left_join(preds, spDf)
library(ggplot2)
ggplot(preds) +
    geom_density(aes(x=sel_prob, col=eventLabel))
preds %>% 
    group_by(eventId, eventLabel, Channel) %>% 
    summarise(avg_sel = mean(sel_prob)) %>% 
    ggplot() +
    # geom_histogram(aes(x=avg_sel, fill=species))
    geom_density(aes(x=avg_sel, col=eventLabel)) + facet_wrap(~Channel)

preds %>% 
    group_by(eventId, species) %>% 
    summarise(avg_sel = mean(sel_prob)) %>% 
    View

# TODO: Anne logging dolphins on PASCAL to finish
# TODO: Anne ran dolphin reject template, will try to tweak thresh to see if
#   that helps on event creation lvel
# TODO: Taiki test different selection removal threshold
# TODO: Taiki adding labels for Ship, Dolphin (inc. species), Pm to "NotBW" events
# TODO: reprocess 576 data - drifts 14, 20

pps <- PAMpalSettings(db=c('../Data/MTC_CV/Databases/PG_2_02_09_CCES_014 - Copy.sqlite3',
                           '../Data/MTC_CV/Databases/PG_2_02_09_CCES_020 - Copy.sqlite3'),
                      binaries = c('../Data/MTC_CV/Binaries/PG_2_02_09_CCES_014/',
                                   '../Data/MTC_CV/Binaries/PG_2_02_09_CCES_020/'),
                      sr_hz=576e3, winLen_sec=.0025, filterfrom_khz=10,
                      filterto_khz=NULL)
data <- processPgDetections(pps, mode='db', id='576fix')
data <- filter(data, Channel == '1')
adWig576 <- processAllWig(data, wl=128, sr=192e3, dir='../Data/MTC_CV/Wiggies288', filter=10, channel=1)
data <- calculateICI(data, time='peakTime')
ici <- getICI(data)
ici <- bind_rows(lapply(ici, function(x) x['All_ici']), .id='eventId')
ici <- rename(ici, ici=All_ici)
adWig576 <- left_join(adWig576, ici)
adWig576$file <- basename(adWig576$file)
adWig576$species <- NULL
write.csv(adWig576, file='../Data/MTC_CV/ADRIFT_Wigner_576fix.csv', row.names=FALSE)

preds576 <- readr::read_csv('../../beaker-wigner-class/preds/adrift_test/adrift_test576/ADRIFT_Wigner_576fix_pred.csv')
table(preds576$db)
table(preds288$db)
View(preds576[!(preds576$UID %in% preds288$UID), ])

(preds288 %>% 
    filter(grepl('CCES_014|CCES_020', db)) %>% 
    ggplot(aes(x=sel_prob, col=eventLabel)) + 
    geom_density()) |
(preds576 %>% 
    ggplot(aes(x=sel_prob, col=eventLabel)) +
    geom_density())

predFixed <- preds288 %>% 
    filter(!grepl('CCES_014|CCES_020', db))
table(basename(predFixed$db))
predFixed <- bind_rows(predFixed, preds576)
write.csv(predFixed, '../../beaker-wigner-class/preds/adrift_test/ADRIFT_Wigner_576Fixed_pred.csv', row.names=FALSE)
rm(preds288, preds576)

# re-run olds to capture missing bins from transfer error and combine
dbs <- list.files('../Data/MTC_CV/Databases/', full.names=TRUE)
dbs <- dbs[c(1:2, 5:10)]
bins <- list.dirs('../data/MTC_CV/Binaries/', full.names=TRUE, recursive = FALSE)
bins <- bins[c(1:2, 5:10)]
pps <- PAMpalSettings(db=dbs, binaries = bins, sr_hz=288e3, filterfrom_khz=10, filterto_khz=NULL, winLen_sec=.0025)
data288 <- processPgDetections(pps, mode='db', id='288again')
data288 <- fitlter(data288, Channel == '1')
data288 <- calculateICI(data288, time='peakTime')
adWig288 <- processAllWig(data288, wl=128, sr=192e3, dir='../Data/MTC_CV/Wiggies288', filter=10, channel=1, retry=FALSE)
sum(grepl('ADRIFT_050.*\\.OE[1-3]{0,1}[0-9]{1}$',names(events(data288))))
drop40 <- grepl('ADRIFT_050.*\\.OE[1-3]{0,1}[0-9]{1}$', names(events(data288)))
drop40 <- drop40 | grepl('ADRIFT_050.*\\.OE40$', names(events(data288)))
names(events(data288))[drop40]
saveRDS(data288, file='../Data/MTC_CV/ADRIFT_Study_288only_withdupes.rds')
data288 <- data288[!drop40]
dataCombined <- bindStudies(data, data288)

table(species(setSpecies(dataCombined, 'pamguard')))
saveRDS(dataCombined, file='../Data/MTC_CV/ADRIFT_Study_CombinedFixed576_DupeDropped.rds')

ad288 <- readr::read_csv('../Data/MTC_CV/ADRIFT_Wigner288.csv')
ad576 <- readr::read_csv('../Data/MTC_CV/ADRIFT_Wigner_576fix.csv')
adComb <- ad288 %>% 
    filter(!grepl('CCES_014|CCES_020', db)) %>% 
    bind_rows(ad576)
write.csv(adComb, file='../Data/MTC_CV/ADRIFT_Wigner_FixedCombined.csv')
