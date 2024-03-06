#### Original New Data Wigner Processing ####
source('../../beaker-wigner-class/R/CV4EProcessingFunctions.R')
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

#### Adding bad noise for selnet retrain ####
getBadSelDets <- function(x, cv, pct=.03, angleCut=90, nMin=1e3, 
                          species = c('ZC', 'NotBW', 'BW43', 'MS', 'BW39V'),
                          seed=112188) {
    if(all(x$angle <= 2*pi)) {
        x$angle <- x$angle * 180 / pi
    }
    x <- group_by(x, eventId) %>% 
        mutate(n=n()) %>% 
        ungroup() %>% 
        filter(n > nMin,
               angle < angleCut,
               eventLabel %in% species)
    cv <- filter(cv, species %in% c('ZC', 'MS', 'BB', 'BW43', 'BW37V'))
    origSp <- table(cv$species)
    newCounts <- round(origSp * pct, 0)
    set.seed(seed)
    bads <- sample_n(x, size=sum(newCounts))
    bads$sm3m <- FALSE
    bads$survey <- 'pascal'
    bads$species <- rep(names(newCounts), newCounts)
    bads$Latitude <- NA
    bads$Longitude <- NA
    select(bads, all_of(c(colnames(cv), 'eventLabel')))
}

cvData <- readr::read_csv('../../beaker-wigner-class/labels/Combined_BW10k_Train.csv')
newData <- readr::read_csv('../Data/MTC_CV/ADRIFT_Wigner_FixedCombined.csv')
bads <- getBadSelDets(newData, cvData, pct=.08, nMin=400)
table(bads$drift, bads$eventLabel)
file.copy(from=file.path('../Data/MTC_CV/Wiggies288', bads$file),
    to=file.path('../../beaker-wigner-class/data/MTCNoise8', bads$file))
bads$file <- paste0('MTCNoise8/', bads$file)
str(bads)
bads$eventLabel <- NULL
withBads <- bind_rows(cvData, bads)
write.csv(withBads, file='../../beaker-wigner-class/labels/CombinedWithBads8_BW10k_Train.csv', row.names = FALSE)
