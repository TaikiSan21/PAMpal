data <- readRDS('../Data/3D2PAMr_test_files/AMStudy.rds')
data <- addRecordings(data, folder='../Data/3D2PAMr_test_files/Recordings/')

clips <- getClipData(data[3], mode='detection')

plot(clips[[1]])


doWig <- function(x, wl=256, dir='../../CV4Ecology/images/') {
    for(i in seq_along(x)) {
        clip <- PAMpal:::clipAroundPeak(x[[i]]@.Data[, 1], length=wl)
        wig <- PAMmisc::wignerTransform(clip, n=wl, sr=x[[i]]@samp.rate, plot=FALSE)
        png(filename = file.path(dir, paste0(names(x)[i], '.png')), width=4, height=4, res=300, units='in')
        par(mar=c(0,0,0,0))
        image(t(wig$tfr), col=viridisLite::viridis(25), useRaster=TRUE)
        dev.off()
    }
}

w <- doWig(clips, wl=128)
oPar <- par(mar=c(0,0,0,0))
image(x=w$t, y=w$f, z=t(w$tfr), col=viridisLite::viridis(25), useRaster=TRUE)
par(oPar)
