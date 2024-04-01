clickSlope <- function(data, calibration=NULL, sFreq=c(10, 20)) {
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
        ##### DO STUFF AFTER HERE #######
        thisWave <- data$wave[, chan]
        
        ####### SAMPLE CALIBRATION CODE, MODIFY AS NEEDED ###
        # Calibration applies to power spec (20*log10 space)
        thisSpec <- spec(thisWave, f=data$sr, wl=512, norm=FALSE, correction='amplitude', plot=FALSE)
        # `spec` has problems with high samplerates producing NA frequencies (integer overflow)
        # so we recreate the frequency vector in a more stable way
        freq <- seq(from=0, by = thisSpec[2,1] - thisSpec[1,1], length.out = nrow(thisSpec))
        dB <- 20*log10(thisSpec[,2])
        roll3 <- data.table::frollmean(dB, n=3)
        roll5 <- data.table::frollmean(dB, n=5)
        # plot(x=freq, y=dB, type='l', col='black')
        # lines(x=freq, y=roll3, col='blue')
        # lines(x=freq, y=roll5, col='purple')
        ix1 <- which.min(abs(freq - sFreq[1]))
        ix2 <- which.min(abs(freq - sFreq[2]))
        slopeDb <- c(dB[ix1], dB[ix2])
        slopeFreq <- c(freq[ix1], freq[ix2])
        slope <- (slopeDb[2] - slopeDb[1]) / (slopeFreq[2] - slopeFreq[1])
        # lines(x=slopeFreq, y=slopeDb, col='red')
        # This chunk searches for the calibration function added with `addCalibration`
        if(!is.null(calibration)) {
            #DO CA
            if(is.function(calibration)) {
                calFun <- calibration
            } else if(is.character(calibration)) {
                calFun <- findCalibration(calibration)
            }
            # Function takes in frequnecy in kHz and returns dB correction
            dB <- dB + calFun(freq * 1e3)
        }
        # dB is now adjusted for calibration, use as you like:
        calibratedMax <- freq[which.max(dB)]
        #### END CALIBRATION SECTION #####
        
        # Store results you want output in `thisResult`:
        # wavMax <- max(thisWave)
        # thisResult$wavMax <- wavMax
        #### STORE THIS CHANNEL ####
        result[[chan]] <- thisResult
    }
    # Combine results for each channel into a single dataframe
    # Theres no need to label each channel, this handles in earlier steps within PAMpal
    result <- bind_rows(result)
    result
}