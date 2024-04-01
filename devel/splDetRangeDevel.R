db <- '../Data/SoundscapeTest/ADRIFT_017.sqlite3'
bin <- '../Data/SoundscapeTest/ADRIFT_017/'
library(PAMpal)
pps <- PAMpalSettings(db, bin, sr_hz='auto', filterfrom_khz=10, filterto_khz=NULL, winLen_sec=.0025)
data <- processPgDetections(pps, id='ADRIFT_017')
saveRDS(data, file = '../Data/SoundscapeTest/adrift17study.rds')

library(readr)
gps <- read_csv('../Data/SoundscapeTest/ADRIFT_017_GPS.csv')
data <- filter(data, UTC >= min(gps$UTC),
               UTC <= max(gps$UTC))
data <- addGps(data, gps=gps)
ol <- read_csv('../Data/SoundscapeTest/metrics/ADRIFT_017_OL_2min.csv')
tol <- read_csv('../Data/SoundscapeTest/metrics/ADRIFT_017_TOL_2min.csv')
tol25 <- read_csv('../Data/SoundscapeTest/metrics/ADRIFT_017_TOL_pct25_2min.csv')
tol25$UTC <- tol25[[1]] + 7 * 3600
names(tol25)[2:24] <- paste0(names(tol25)[2:24], '_25p')
bb <- read_csv('../Data/SoundscapeTest/metrics/ADRIFT_017_BB_2min.csv')
# these were not originally in UTC
bb$UTC <- bb[[1]] + 7 * 3600
# bb[[1]] <- NULL

data <- matchTimeData(data, bb[c(2,3)], mode='event', thresh=3600, replace=TRUE)
data <- matchTimeData(data, tol25[c('UTC', 'TOL_1000_25p')], mode='event', thresh=3600, replace=TRUE)
data <- matchEnvData(data, nc='etopo180', var='altitude')
hm <- getMeasures(data)
hm$nClicks <- sapply(events(data), nClicks)
plot(hm$`BB_100-24000`)
library(ggplot2)

ggplot(hm) +
    geom_point(aes(x=`BB_100-24000`, y=log10(nClicks+1)))

ggplot(hm) +
    geom_point(aes(x=abs(altitude_mean), y=`BB_100-24000`))


estimateRange = function(f=500, h=50, SL=120, NL=60, r=1, SNRthresh=1){
    
    # freq dep absorption
    alpha = AcousticAbsorption(f)
    SNRcalc=500
    while(SNRcalc>SNRthresh){
        SNRcalc = SL - (10log10(r)+10log10(h/2)+(alpha/1000)*r) - NL
        r = r+1
    }
    return(r)
}

alpha <- .06 #for 1khz .02 for 500
logfun <- function(SL=120, NL=60, SNR=1, h=50, alpha=.02) {
    function(x) {
        abs(SL - 10 * log10(x) -10*log10(h/2) - x/1000*alpha - NL - SNR)
    }
}
optim(par=10e3, fn=logfun(120,70,10), method='Brent', lower=1, upper=100e3)

hm$estRange <- NA
for(i in 1:nrow(hm)) {
    if(is.na(hm$TOL_1000_25p[i])) {
        next
    }
    opTry <- optim(par=10e3, fn=logfun(150, hm$TOL_1000_25p[i], 1, abs(hm$altitude_mean[i]), .02),
                   method='Brent', lower=1, upper=500e3)
    hm$estRange[i] <- opTry$par
}
plot(hm$estRange)
hm$ix <- 1:nrow(hm)
hm$nClicks <- sapply(events(data), nClicks)

library(patchwork)
(ggplot(hm) +
    geom_point(aes(x=ix, y=estRange)) + 
        ggplot(hm) +
        geom_point(aes(x=ix, y=log10(nClicks)))) /
(ggplot(hm, aes(x=ix)) +
    geom_point(aes(y=altitude_mean)) +
        ggplot(hm) +
        geom_point(aes(x=ix, y=TOL_1000_25p)))
    
