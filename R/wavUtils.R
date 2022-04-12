# utils for wav file stuff

# get spectrum or cepstrum of a signal
myGram <- function(x, channel=1, wl = 512, window = TRUE, sr=NULL,
                   from=0, to=NULL, noise=FALSE, mode=c('spec', 'ceps'),
                   decimate=1) {
    if(is.list(x)) {
        wave <- x$wave[, channel]
        if(is.null(sr)) {
            sr <- x$sr
        }
    } else if(is.vector(x)) {
        wave <- x
    } else if(inherits(x, 'Wave')) {
        wave <- switch(channel,
                       '1' = x@left,
                       '2' = x@right,
                       return(NULL)
        )
        if(is.null(sr)) {
            sr <- x@samp.rate
        }
    } else if(inherits(x, 'WaveMC')) {
        if(channel > ncol(x@.Data)) {
            return(NULL)
        }
        wave <- x@.Data[, channel]
        if(is.null(sr)) {
            sr <- x@samp.rate
        }
    }
    if(decimate > 1) {
        wave <- decimate(wave, q=decimate)
        sr <- sr / decimate
    }
    if(from > 0) {
        # kinda janky because NULL * 1e3 is not NULL anymore, its numeric(0)
        if(!is.null(to)) {
            to_hz <- to * 1e3
        } else {
            to_hz <- NULL
        }
        wave <- bwfilter(wave, f=sr, n=4, from=from*1e3, to = to_hz, output='sample')
    }
    wave <- clipAroundPeak(wave, wl, noise=noise)
    if(window) {
        wave <- wave * hanning(length(wave)) / mean(hanning(length(wave)))
    }

    switch(match.arg(mode),
           'spec' = {
               FUN <- function(x) {
                   result <- Mod(fft(x))
                   20*log10(result)
               }
               y <- (1:(wl)) / wl * sr
           },
           'ceps' = {
               FUN <- function(x) {
                   result <- Mod(fft(x))^2
                   result <- ifelse(result == 0, 1e-15, result)
                   result <- log(result)
                   # browser()
                   result <- result - mean(result)
                   result <- fft(result, inverse=TRUE)
                   # abs(Re(result))

                   abs(Re(result))
               }
               y <- (1:wl) / sr
           })

    ans <- matrix(NA, nrow=wl%/%2, ncol=2)
    ans[, 1] <- y[1:(wl%/%2)]
    dB <- FUN(wave)[1:(wl)]
    isInf <- is.infinite(dB)
    if(all(isInf)) {
        return(ans)
    }
    if(any(isInf)) {
        dB[dB == Inf] <- max(dB[!isInf], na.rm=TRUE)
        dB[dB == -Inf] <- min(dB[!isInf], na.rm=TRUE)
    }
    ans[, 2] <- dB[1:(wl%/%2)]
    # list(dB = dB[1:(wl%/%2)], freq=y[1:(wl%/%2)])
    ans
}

# list of binary files to matrix of spectra
binToSpecMat <- function(bin, channel=1, freq, wl, filterfrom_khz, filterto_khz,
                         sr, decimate, mode, noise=FALSE, calFun, ...) {
    lapply(bin, function(x) {
        if(is.null(x$wave)) {
            return(rep(NA, wl%/%2))
        }
        hasChan <- 1:ncol(x$wave)
        useChan <- hasChan[hasChan %in% channel]
        chanTemp <- matrix(NA, nrow=length(freq), ncol=length(useChan))
        for(c in seq_along(useChan)) {
            chanTemp[, c] <- myGram(x, channel = useChan[c], wl=wl,
                                    from=filterfrom_khz,
                                    to=filterto_khz,
                                    sr = sr, noise=noise,
                                    decimate = decimate,
                                    mode=mode, ...)[, 2] + calFun(freq)
        }
        chanTemp <- chanTemp[, apply(chanTemp, 2, function(x) !any(is.na(x)))]
        if(all(chanTemp == 0)) {
            return(rep(NA, wl%/%2))
        }
        doLogAvg(chanTemp, log = mode=='spec')
    })
}

# these two adjus 0-255 scale matrices for brightness/contrast
doBrightness <- function(x, val=0) {
    if(val == 0) {
        return(x)
    }
    val <- val * -1
    x <- x + val
    x[x > 255] <- 255
    x[x < 0] <- 0
    x
}

doContrast <- function(x, val=0) {
    if(val == 0) {
        return(x)
    }
    factor <- (259 * (val + 255)) / (255 * (259 - val))
    x <- factor * (x - 128) + 128
    x[x > 255] <- 255
    x[x < 0] <- 0
    x
}

# Wave class object or vector waveform to a spectrogram of signal
wavToGram <- function(wav, sr=NULL, wl=1024, hop=.5, mode=c('spec', 'ceps'), axes=FALSE, ...) {
    if(inherits(wav, 'Wave')) {
        if(is.null(sr)) {
            sr <- wav@samp.rate
        }
        wav <- wav@left
    } else if(inherits(wav, 'WaveMC')) {
        if(is.null(sr)) {
            sr <- wav@samp.rate
        }
        wav <- wav@.Data[, 1]
    }
    mode <- match.arg(mode)
    if(hop <= 1) {
        hop <- wl * hop
    }
    nSlices <- ceiling((length(wav) - wl)/(hop)) + 1
    slices <- seq(from=1, by=hop, length.out=nSlices)
    mat <- t(apply(as.matrix(slices), 1, function(s) {
        thisWav <- wav[s:(s+wl-1)]
        thisGram <- myGram(thisWav, mode=mode, wl=wl, sr=sr, ...)
        thisGram[, 2]
    }))
    if(!axes) {
        return(mat)
    }
    yAxis <- myGram(wav[1:wl], mode=mode, wl=wl, sr=sr, ...)[, 1]
    switch(mode,
           'spec' = {
               yAxis <- yAxis / 1e3
           },
           'ceps' = {
               yAxis <- yAxis * 1e3
           }
    )
    xAxis <- slices / sr
    list(mat=mat, x=xAxis, y=yAxis)
}

# clip of fixed length, zeropads if needed and deals with edge case
clipAroundPeak <- function(wave, length, noise=FALSE) {
    if(length(wave) < length) {
        return(c(wave, rep(0, length - length(wave))))
    }
    peak <- which.max(abs(wave))
    low <- peak - floor(length/2)
    high <- ceiling(peak + length/2) -1
    if(low < 1) {
        return(wave[1:length])
    }
    if(high > length(wave)) {
        return(wave[(length(wave)-length+1):length(wave)])
    }
    if(noise) {
        if(low - length > 1) {
            return(wave[(low-length):(low-1)])
        } else {
            return(wave[1:length])
        }
    }
    wave[low:high]
}
