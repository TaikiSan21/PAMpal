

# hist(ici$ici[ici$ici < 2], breaks=50)


doFancyIci <- function(x, time='peakTime', maxICI=2, res=.001, mode=c('ac', 'ceps', 'acslow'), plot=FALSE) {
    x$timeNum <- as.numeric(x[[time]])
    x$timeNum <- round((x$timeNum - min(x$timeNum))/res, 0) * res
    # browser()
    timeSeq <- seq(from=0, to=max(x$timeNum), by=res)
    detSeq <- rep(0, length(timeSeq))
    detSeq[timeSeq %in% x$timeNum] <- 1
    fftLen <- 2^(ceiling(log2(maxICI * 2/res)))
    winSeq <- seq(from=1, to=length(detSeq), by= fftLen / 2)
    outMat <- matrix(0, nrow=fftLen, ncol=length(winSeq))
    allZero <- numeric()
    last <- max(winSeq) + fftLen - 1
    detSeq <- c(detSeq, rep(0, last - length(detSeq)))
    FUN <- switch(match.arg(mode),
                     'ac' = {
                         function(s) {
                             s <- s - mean(s)
                             res <- fft(s)
                             res <- Mod(res)
                             res <- res - mean(res)
                             Re(fft(res, inverse=TRUE))
                         }
                     },
                     'ceps' = {
                         function(s) {
                             res <- Mod(fft(s))^2
                             res <- ifelse(res == 0, 1e-15, res)
                             res <- log(res)
                             res <- res - mean(res)
                             res <- fft(res, inverse=TRUE)
                             # abs(Re(res))
                             Re(res)
                         }
                     },
                     'acslow' = {
                         function(s) s
                     }
    )
    for(i in seq_along(winSeq)) {
        thisWin <- detSeq[winSeq[i]:(winSeq[i]+fftLen-1)]
        if(sum(thisWin) < 2) {
            allZero <- c(allZero, i)
            next
        }
        outMat[, i] <- FUN(thisWin)
    }
    outMat <- outMat[, -allZero]
    avg <- apply(outMat, 1, mean)[1:(fftLen/2)]
    max <- which.max(avg)
    if(plot) {
        plot(x=timeSeq[1:(fftLen/2)], y=avg, type='l')
        title(main=paste0('Mode = ', mode, '  Max = ', round(timeSeq[max], res^-1)))
    }
    list(mat=outMat, avg=avg, max=timeSeq[max])
}

data <- readRDS('../Data/3D2PAMr_test_files/AMStudy.rds')
library(PAMpal)
data <- calculateICI(data, time = 'peakTime')
ici <- getICI(data[[4]], 'data')
ici <- ici[ici$detectorName == 'Click_Detector_1',]

hm <- doFancyIci(ici, maxICI=2, res=.001, plot=T, mode='ceps')
hm <- doFancyIci(ici, maxICI=.25, res=.001, plot=T, mode='ac')

fake <- c(seq(from=1.03, length.out=20, by=.674),
          seq(from=4.52, length.out=15, by=.41))
fake <- data.frame(peakTime=fake)
fake$peakTime <- fake$peakTime + rnorm(nrow(fake), mean=0, sd=.01)
fake$peakTime <- sort(fake$peakTime)
maxICI <- 1
hm <- doFancyIci(fake, maxICI=maxICI, res=.001, plot=T, mode='ac')
hm <- doFancyIci(fake, maxICI=maxICI, res=.001, plot=T, mode='ceps')
hist(diff(fake$peakTime), breaks=30, xlim=c(0, maxICI))
