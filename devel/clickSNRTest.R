# run cces note:
# station 4 has diff duty cyclyes - 2/2 for 20 days, 2/18 for rest
# can deal with this by having separate station IDs for first 20 and rest
# need to manually adjust

# can convert deployment details worksheet to rec info:
dd <- read.csv('../Data/RoboJ/Test/deployDetails.csv', stringsAsFactors = FALSE)

ri <- select(dd, all_of(c('Project', 'Station'='DeploymentID', 'Recorder'='Type', 'deployTime'='Data_Start', 'retrTime'='Data_End',
             'dutyOn'='RecordingDuration_m', 'dutyOff'='RecordingInterval_m')))




sqrt(sum(x^2)) / sqrt(sum(y^2))

powsig <- sum(x^2)/n
pownoise <- sum(y^2)/n

snr <- powsig/pownoise
snrdb <- 10*log10(snr)
click <- loadPamguardBinaryFile('../Data/3D2PAMr_test_files/BinaryFiles_20160815_WithUID/Click_Detector_Click_Detector_Clicks_20160815_200001.pgdf')
click <- loadPamguardBinaryFile('../Data/RoboJ/Test/Binaries_Station-7_Card-A_WithUID/20160824/Click_Detector_Click_Detector_Clicks_20160824_000207.pgdf')
click <- click$data

doSnr <- function(x, wSig=256, wNoise=NULL) {
    if(is.null(wNoise)) {
        wNoise <- wSig
    }
    if(is.list(x) &&
       'wave' %in% names(x)) {
        x <- x$wave[, 1]
    }
    sig <- PAMpal:::clipAroundPeak(x, wSig, FALSE)
    noise <- PAMpal:::clipAroundPeak(x, wNoise, TRUE)
    # not sqrt because we are just logging below
    sigRMS <- mean(sig^2)
    noiseRMS <- mean(noise^2)
    10*log10(sigRMS/noiseRMS)
}
wSig <- 256
wNoise <- 1024
doSnr(click[[1]], wSig, wNoise)
sapply(click[1:10], doSnr)

par(mfrow=c(3, 4))
for(i in 1:12) {
    wv <- click[[i+12*2]]$wave[, 2]
    # wv <- PAMpal:::clipAroundPeak(wv, wl)
    sig <- PAMpal:::clipAroundPeak(wv, wSig, FALSE, TRUE)
    # sig <- doIntEnergy(wv, pct=.9)
    noise <- PAMpal:::clipAroundPeak(wv, wNoise, TRUE, TRUE)
    plot(wv, type='l')
    lines(x=sig$ix, y=sig$wav, col='blue')
    snr <- round(doSnr(wv, wSig=wSig, wNoise=wNoise), 2)
    # lines(x=rep(wl, 2), y=c(-1,1), col='blue', lty=3)
    lines(x=noise$ix, y=noise$wav, col='red')
    title(paste0('WL ', wSig, ', SNR ', snr))
}

doIntEnergy <- function(wav, pct=.9) {
    intE <- cumsum(wav^2)
    total <- tail(intE, n=1)
    # q <- (1-pct)/2
    begin <- tail(which(intE < (.5 * total - pct*total/2)), n=1)
    finish <- tail(which(intE < (.5 * total + pct*total/2)), n=1)
    list(ix=begin:finish, wav=wav[begin:finish])
}
intE<-cumsum(wv^2)
totalE<-tail(intE,n=1)
begin<-tail(which(intE<(.5*totalE - (.95*totalE/2))), n=1)
finish<-tail(which(intE<(.5*totalE + (.95*totalE/2))), n=1)
dur<-(finish - begin)/sr

