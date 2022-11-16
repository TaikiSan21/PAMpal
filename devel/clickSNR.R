clickSNR <- function(data, wSig=256, wNoise=NULL, calibration=NULL) {
    result <- list()
    if(is.null(wNoise)) {
        wNoise <- wSig
    }
    # Do same stuff for each channel
    for(chan in 1:ncol(data$wave)) {
        # create storage for this channels outputs
        thisResult <- list()
        ##### DO STUFF AFTER HERE #######
        thisWave <- data$wave[, chan]
        sig <- PAMpal:::clipAroundPeak(thisWave, wSig, FALSE)
        noise <- PAMpal:::clipAroundPeak(thisWave, wNoise, TRUE)
        # not sqrt because we are just logging below
        thisResult$sigRMS <- 10*log10(mean(sig^2))
        thisResult$noiseRMS <- 10*log10(mean(noise^2))
        thisResult$SNR <- thisResult$sigRMS - thisResult$noiseRMS
        #### STORE THIS CHANNEL ####
        result[[chan]] <- thisResult
    }
    # Combine results for each channel into a single dataframe
    # Theres no need to label each channel, this handles in earlier steps within PAMpal
    result <- bind_rows(result)
    result
}
##### Example how to add function, just need to source first or you'll get "cannot find function clickSNR" ####
# db <- '../Data/RoboJ/Test/Station-7_Card-A_MASTER-BW - Copy.sqlite3'
# bin <- '../Data/RoboJ/Test/Binaries_Station-7_Card-A_WithUID/'
# pps <- PAMpalSettings(db,bin)
# set wSig and wNoise to whatever you want - window lengths for measuring sig/noise levels
# pps <- addFunction(pps, clickSNR, module='ClickDetector', wSig=256, wNoise=1024)
#
# data <- processPgDetections(pps, mode='db')
# Now you'll have sigRMS, noiseRMS, and SNR in outputs
