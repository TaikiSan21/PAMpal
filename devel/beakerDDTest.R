# process bw data and testing SSP raytracing for issues

bwdb <- '../Data/3D2PAMr_test_files/HB1603_MF_Master_Leg3_Skala.sqlite3'
bwbin <- '../Data/3D2PAMr_test_files/BinaryFiles_20160815_WithUID/'
bwpps <- PAMpalSettings(bwdb, bwbin, sr_hz=192e3, filterfrom_khz=10, filterto_khz=NULL, winLen_sec=.0025)
bwdata <- processPgDetections(bwpps)
bwdata <- bwdata %>% 
    addGps() %>% 
    addHydrophoneDepth(depthCol='Sensor_0_Depth')


wavDir <- '../Data/3D2PAMr_test_files/WavClips/'
wavList <- list.files(wavDir, pattern='wav$', full.names=TRUE)
bwParams <- list(freqLow=25e3, 
                 freqHigh=70e3,
                 minDepth=500,
                 maxRange=4e3)
bwdata <- calcDiveDepth(bwdata, 
                        wav=wavDir,
                        params=bwParams, 
                        nPlot=20,
                        hpDepthError = 3,
                        outDir='../Data/3D2PAMr_test_files/Outputs')
wat <- filterClicks(bwdata, time=10, speed=50/30)

ggplot(wat, aes(x=UTC, y=-maxDepth, col=keepClick)) +
    geom_point()
maybe <- runDepthReview(wat)



bwclicks <- getClickData(bwdata)

bwclicks <- dplyr::filter(bwclicks, grepl('117', eventId))

depths <- matchEnvData(bwclicks[1,], nc='HYCOM', raw=TRUE, timeout=360, depth=0:500)
library(PAMmisc)
ssp <- createSSP(bwclicks[1, ])
gamout <- mgcv::gam(formula=speed ~ s(depth, k=min(30, length(ssp[[1]]$speed))), data=ssp[[1]])
newdata <- data.frame(depth= 0:max(ssp[[1]]$depth))
predicted <- as.numeric(predict.gam(gamout, newdata=newdata))
newdata$soundSpeed <- predicted
plot(x=ssp[[1]]$speed, y=-ssp[[1]]$depth, type='l', ylim=c(-50, 0))

rt <- raytrace(1000, 1000, c(-85,-85.5),  1, zz=newdata$depth, cc=newdata$soundSpeed, plot=T)
sapply(rt$z, function(x) {
    min(which(x == 0))
})
diffs <- sapply(length(rt$x[[1]]), function(x) {
    sqrt( (rt$x[[1]][x]-rt$x[[2]])^2 + (rt$z[[1]][x] - rt$z[[2]][x])^2)
})
d <- 1000
z <- 8
x <- 2000
(atan((d+z)/x) - atan((d-z)/x))*180/pi
