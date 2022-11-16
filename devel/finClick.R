
finClick <- function(data, sr_hz='auto', winLen_sec=1, threshold=-5, calibration=NULL) {
    # List names of all required packages here to guarantee they will be
    # installed when others try to use your function, ex:
    # packageList <- c('seewave', 'signal')
    packageList <- c()
    for(p in packageList) {
        if(!require(p, character.only = TRUE)) {
            install.packages(p)
            require(p, character.only = TRUE)
        }
    }
    # prepping output
    result <- list()
    # Do same stuff for each channel
    for(chan in 1:ncol(data$wave)) {
        # create storage for this channels outputs
        thisResult <- list()
        if(sr_hz == 'auto') {
            sr <- data$sr
        } else {
            sr <- sr_hz
        }
        # browser()
        ##### DO STUFF AFTER HERE #######
        thisWave <- data$wave[, chan]
        fftSize <- round(sr * winLen_sec, 0)
        fftSize <- fftSize + (fftSize %% 2)
        wherePeak <- which.max(abs(thisWave))
        starts <- seq(from=-.5, to=.5, by=.05) * fftSize + wherePeak - fftSize/4
        windows <- vector('list', length=length(starts))
        for(s in seq_along(starts)) {
            start <- max(starts[s], 1)
            end <- min(starts[s] + fftSize / 2, length(thisWave))
            if(start > length(thisWave) ||
               end < 1) {
                windows[[s]] <- list(dB=NA_real_, freq=0, time=starts[s]/sr)
                next
            }
            fftVec <- start:end
            winWave <- thisWave[start:end]
            if(length(winWave) < fftSize/2) {
                winWave <- c(winWave, rep(0, fftSize/2 - length(winWave)))
            }
            if(all(winWave == 0)) {
                windows[[s]] <- list(dB=NA_real_, freq=0, time=starts[s]/sr)
                next
            }
            spec <- seewave::spec(winWave, wl=fftSize/2, f=sr, norm=FALSE, correction='amplitude', plot=FALSE)
            dB <- 20*log10(spec[, 2])
            whichMax <- which.max(dB)[1]
            windows[[s]] <- list(dB=dB[whichMax], freq=spec[whichMax, 1], time=starts[s]/sr)
        }
        windows <- bind_rows(windows)
        windows$dB[is.na(windows$dB)] <- min(windows$dB, na.rm=TRUE) - 10
        windows$dB <- windows$dB - max(windows$dB)
        useCalc <- windows[windows$dB > threshold, ]
        peak <- useCalc$freq[which.max(useCalc$dB)[1]]
        center <- mean(range(useCalc$freq))
        if(nrow(useCalc) > 1) {
            lm <- stats::lm(freq ~ time, data=useCalc)
            slope <- unname(lm$coefficients[2])
        } else {
            slope <- 0
        }
        # thisSpec <- spec(thisWave, f=data$sr, wl=512, norm=FALSE, correction='amplitude', plot=FALSE)
        # `spec` has problems with high samplerates producing NA frequencies (integer overflow)
        # so we recreate the frequency vector in a more stable way
        # freq <- seq(from=0, by = thisSpec[2,1] - thisSpec[1,1], length.out = nrow(thisSpec))
        # dB <- 20*log10(thisSpec[,2])

        # Store results you want output in `thisResult`:
        # wavMax <- max(thisWave)
        # thisResult$wavMax <- wavMax
        #### STORE THIS CHANNEL ####
        result[[chan]] <- list(freqSlope=slope, centerFreq=center, cenPeakDiff=unname(center-peak))
    }
    # Combine results for each channel into a single dataframe
    # Theres no need to label each channel, this handles in earlier steps within PAMpal
    result <- bind_rows(result)
    result
    # starts
}
