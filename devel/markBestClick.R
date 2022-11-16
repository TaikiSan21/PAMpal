markBestClick <- function(data, calibration=NULL, sr_hz='auto', winLen_sec=.0025, from_khz=0, to_khz=NULL) {
  packageList <- c('seewave')
  for(p in packageList) {
    if(!require(p, character.only = TRUE)) {
      install.packages(p)
      require(p, character.only = TRUE)
    }
  }
  if(sr_hz == 'auto') {
    sr <- data$sr
  } else {
    sr <- sr_hz
  }
  fftLen <- round(winLen_sec * sr, 0)
  
  result <- list()
  for(chan in 1:ncol(data$wave)) {
    thisResult <- list(bestClick=FALSE)
    thisWave <- data$wave[, chan]
    thisWave <- PAMpal:::clipAroundPeak(thisWave, length = fftLen)
    
    ####### SAMPLE CALIBRATION CODE, MODIFY AS NEEDED ###
    # Calibration applies to power spec (20*log10 space)
    thisSpec <- spec(thisWave, f=data$sr, wl=fftLen, norm=FALSE, correction='amplitude', plot=FALSE)
    # `spec` has problems with high samplerates producing NA frequencies (integer overflow)
    # so we recreate the frequency vector in a more stable way
    thisSpec[, 1] <- seq(from=0, by = thisSpec[2,1] - thisSpec[1,1], length.out = nrow(thisSpec))
    thisSpec[, 2] <- 20*log10(thisSpec[,2])
    # This chunk searches for the calibration function added with `addCalibration`
    if(!is.null(calibration)) {
      #DO CA
      if(is.function(calibration)) {
        calFun <- calibration
      } else if(is.character(calibration)) {
        calFun <- findCalibration(calibration)
      }
      # Function takes in frequnecy in kHz and returns dB correction
      thisSpec[, 2] <- thisSpec[, 2] + calFun(freq * 1e3)
    }
    aboveLow <- thisSpec[, 1] >= from_khz
    if(is.null(to_khz)) {
      underHigh <- rep(TRUE, nrow(thisSpec))
    } else {
      underHigh <- thisSpec[, 1] <= to_khz
    }
    inRange <- thisSpec[aboveLow & underHigh, ]
    thisResult$maxdB <- max(inRange[, 2])
    thisResult$maxFreq <- inRange[which.max(inRange[, 2]), 1]
    result[[chan]] <- thisResult
  }
  result <- bind_rows(result)
  result$bestClick[which.max(result$maxdB)] <- TRUE
  result$maxFreq <- NULL
  result$maxdB <- NULL
  result
}
