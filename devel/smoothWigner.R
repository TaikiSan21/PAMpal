smoothWVD <- function(signal, n=NULL, sr, plot=FALSE, tw=NULL, fw=NULL) {
    if(inherits(signal, 'Wave')) {
        sr <- signal@samp.rate
        signal <- signal@left / 2^(signal@bit - 1)
    }
    if(inherits(signal, 'WaveMC')) {
        sr <- signal@samp.rate
        signal <- signal@.Data[, 1] / 2^(signal@bit - 1)
    }
    if(is.null(tw)) {
        tw <- n %/% 10
    }
    tw <- tw + 1 - (tw %% 2)
    tw <- signal::hamming(tw)
    if(is.null(fw)) {
        fw <- n %/% 4
    }
    fw <- fw + 1 - (fw %% 2)
    fw <- signal::hamming(fw)
    tmid <- (length(tw)-1) %/% 2
    fmid <- (length(fw)-1) %/% 2

    analytic <- PAMmisc:::toAnalytic(signal)[1:length(signal)] # size changed during toAnalytic function
    conjAnalytic <- Conj(analytic)
    if(is.null(n)) {
        n <- PAMmisc:::nextExp2(length(analytic))
    }

    nRow <- n # nFreq bins
    nCol <- length(analytic) # nTimesteps
    # nCol <- n # nTimesteps

    tfr <- matrix(0, nRow, nCol)

    for(iCol in 1:nCol) {
        taumax <- min(iCol-1, nCol-iCol, round(nRow/2)-1, fmid)
        # cat('Max', iCol + taumax,
        #     'Min', iCol-taumax, '\n')
        tau <- -taumax:taumax
        indices <- (nRow + tau) %% nRow + 1
        # * .5 in PG?
        # print(tau+fmid)
        # browser()
        tfr[indices, iCol] <- fw[tau + fmid+1] * analytic[iCol+tau] * conjAnalytic[iCol-tau] / 2

        tau <- round(nRow/2)
        if(iCol + tau <= nCol &&
           iCol - tau >= 1) {
            # browser()
            # PG is like this, wv.wge is just the same fucking thing???
            tfr[tau+1, iCol] <- fw[tau+fmid+1]*(analytic[iCol+tau] * conjAnalytic[iCol-tau] +
                                     analytic[iCol-tau] * conjAnalytic[iCol+tau])/2
        }
    }
    tfr <- apply(tfr, 2, fft)
    result <- list(tfr=Re(tfr), t=1:nCol/sr, f=sr/2*1:nRow/nRow)
    if(plot) {
        image(t(result$tfr), xaxt='n', yaxt='n',
              ylab='Frequency (kHz)', xlab = 'Time (ms)',
              col = viridisLite::viridis(25), useRaster=TRUE)
        xPretty <- pretty(result$t, n=5)
        # axis(1, at = 1:4/4, labels = round(1e3*max(result$t)*1:4/4, 3))
        axis(1, at=xPretty / max(result$t), labels=xPretty*1e3)
        yPretty <- pretty(result$f, n=5)
        # axis(2, at = 1:4/4, labels = round(max(result$f)*1:4/4/1e3, 1))
        axis(2, at = yPretty / max(result$f), labels=yPretty/1e3)
    }
    result
}

swig <- smoothWVD(clip, n=128, sr=192e3, plot=T, fw=33)
wig <- PAMmisc::wignerTransform(clip, n=128, sr=192e3, plot=T)

bin <- getBinaryData(bwData, 2082000121)
clip <- PAMpal:::clipAroundPeak(bin[[1]]$wave[, 1], 128)
