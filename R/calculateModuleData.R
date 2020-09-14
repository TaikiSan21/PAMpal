#' @title Run Custom Calculations on Pamguard Module Data
#'
#' @description Run a list of custom calculations on a Pamguard binary
#'   file.
#'
#' @param binData Pamguard binary data as read in by
#'   \code{\link[PamBinaries]{loadPamguardBinaryFile}}
#' @param binFuns A named list of functions to run on each Pamguard module.
#'   Currently supported modules are 'ClickDetector' and 'WhistlesMoans', a
#'   sample input for binFuns would be list('ClickDetector'=list(clickFun1,
#'   clickFun2), 'WhistlesMoans'=list(wmFun1))
#'
#' @return A data frame with one row for each channel of each detection.
#'   Each row will have the UID, channel number, and name of the detector.
#'   Clicks of different classifications are treated as different detectors
#'   for this purpose, with the classification label number appended to the
#'   detector name. The number of columns will depend on the results of the
#'   calculations from the supplied binFuns.
#'
#' @author Taiki Sakai \email{taiki.sakai@@noaa.gov}
#'
#' @import dplyr
#' @importFrom PamBinaries convertPgDate contourToFreq
#' @export
#'
calculateModuleData <- function(binData, binFuns=list('ClickDetector'=list(standardClickCalcs))) {
    if(length(binData$data) == 0 ||
       is.null(binData$data[[1]]$UID)) {
        return(NULL)
    }
    moduleType <- binData$fileInfo$fileHeader$moduleType
    detName <- binData$fileInfo$fileHeader$moduleName
    moduleType <- gsub(' ', '', moduleType)
    # TEMP FIX FOR NOW, MAYBE IS FINE?
    if(moduleType == 'SoundTrapClickDetector') {
        moduleType <- 'ClickDetector'
    }
    if(moduleType == 'ClickDetector' &&
       binData$fileInfo$fileHeader$streamName == 'Trigger Background') {
        moduleType <- 'ClickDetectorTriggerBackground'
        detName <- paste0(detName, '_TriggerBackground')
    }
    # For now, cepstrum-ing like dis
    if(moduleType == 'WhistlesMoans' &&
       grepl('cepstrum', detName, ignore.case = TRUE)) {
        moduleType <- 'Cepstrum'
    }
    detName <- gsub(' ', '_', detName)
    if(!(moduleType %in% names(binFuns)) ||
       length(binFuns[[moduleType]])==0) {
        warning("I don't have functions for Module Type ", moduleType)
        # If nothing, just UID and detectorName? Fine for now
        result <- data.frame(UID = sapply(binData$data, function(x) x$UID),
                             UTC = sapply(binData$data, function(x) x$date))
        result$UID <- as.character(result$UID)
        result$UTC <- convertPgDate(result$UTC)
        result$detectorName <- rep(detName, nrow(result))
        result$BinaryFile <- rep(basename(binData$fileInfo$fileName), nrow(result))
        return(result)
    }
    # Adding this to binFuns for PG version so we always at least get UID and time
    # ClickDetector gets 1 row for each channel, has 'nChan' in it
    getBasic <- function(type) {
        switch(type,
               'ClickDetector' = function(x) {
                   nRep <- if('nChan' %in% names(x)) {
                       x$nChan
                   } else {
                       1
                   }
                   result <- data.frame(UID = rep(as.character(x$UID), nRep),
                                        UTC = rep(x$date, nRep),
                                        stringsAsFactors = FALSE)
                   angs <- x$angles
                   if(length(angs) == 2) {
                       warning('More than one angle for click UID' , result$UID[1],
                               ', only using first.')
                   }
                   result$angle <- angs[1]
                   result$angleError <- x$angleErrors[1]
                   result
               }, function(x) {
                   nRep <- if('nChan' %in% names(x)) {
                       x$nChan
                   } else {
                       1
                   }
                   data.frame(UID = rep(as.character(x$UID), nRep),
                              UTC = rep(x$date, nRep),
                              stringsAsFactors = FALSE)
               })

    }
    # if you add new modules need to add to addBinaries - weve filtered them out there
    result <- switch(
        moduleType,
        'ClickDetector' = {
            allClicks <- doClickCalcs(binData$data, c(getBasic(moduleType), binFuns[['ClickDetector']]))
            # We want each 'type' of click to be separate 'detector'. This is PG only.
            allNames <- bind_rows(
                lapply(binData$data[unique(as.character(allClicks$UID))], function(x) {
                    data.frame(UID=as.character(x$UID),
                               detectorName=unique(c(x$type, unlist(x$annotations$classification))),
                               stringsAsFactors = FALSE)
                })) %>%
                mutate(detectorName = paste(detName, .data$detectorName, sep='_'))
            allClicks <- left_join(allClicks, allNames, by='UID')
            allClicks$callType <- 'click'
            allClicks
        },
        'WhistlesMoans' = {
            allWhistles <- doWhistleCalcs(binData$data, c(getBasic(moduleType), binFuns[['WhistlesMoans']]))
            allWhistles$detectorName <- detName
            allWhistles$callType <- 'whistle'
            allWhistles
        },
        'Cepstrum' = {
            allCepstrum <- doCepstrumCalcs(binData$data, c(getBasic(moduleType), binFuns[['Cepstrum']]))
            allCepstrum$detectorName <- detName
            allCepstrum$callType <- 'cepstrum'
            allCepstrum
        },
        warning("I don't know how to deal with Module Type ", moduleType)
    )
    result$BinaryFile <- basename(binData$fileInfo$fileName)
    result$UTC <- convertPgDate(result$UTC)
    result
}

# In general just make sure data has $wave and $sr for clicks?
doClickCalcs <- function(clickData, clickFuns) {
    allClicks <- vector('list', length = length(clickFuns))
    # This just returns a df we will bind to db by UID
    for(f in seq_along(clickFuns)) {
        # Apply this function to each datapoint in binary
        allClicks[[f]] <- bind_rows(
            lapply(clickData, function(oneClick) {
                tryCatch({
                    clickFuns[[f]](oneClick)
                }, error = function(e) {
                    # browser()
                    cat('Error in function ', names(clickFuns)[f],
                        ': ', as.character(e$call[1]), '-', e$message,
                        '. UID: ', oneClick$UID, sep='')
                }
                )
            })
        )
    }
    bind_cols(allClicks)
}


# WHISTLES DAKINE
# Hop = (StartSample + 1) / (1st Slice Num)
# FFT Length = Samp Dur - (N slices - 1) * Hop(above)
# FFT Length = SampleRate * max(peakData) / MaxFreq
# Mult bin by SR / Len

# PROBLEM right now this is very PG specific. Which is maybe okay, because
# moduleCalcs is PG-based? Possbly rename these doPgWhistleCalcs or w/e
doWhistleCalcs <- function(whistleData, whistleFuns) {
    # Does this die if 1st slice is #1 (or 0, w/e index is)
    # probably, so use next whistle
    whistleData <- contourToFreq(whistleData)
    # tempData <- whistleData[[1]]
    # if(tempData$sliceData[[1]]$sliceNumber == 0) {
    #     tempData <- whistleData[[2]]
    # }
    # fftHop <- (tempData$startSample + 1)/tempData$sliceData[[1]]$sliceNumber
    # fftLen <- tempData$sampleDuration -
    #     (tempData$sliceData[[tempData$nSlices]]$sliceNumber - tempData$sliceData[[1]]$sliceNumber) * fftHop

    allWhistles <- vector('list', length=length(whistleFuns))
    for(f in seq_along(whistleFuns)) {
        allWhistles[[f]] <- bind_rows(
            lapply(whistleData, function(oneWhistle) {
                tryCatch({
                    whistleFuns[[f]](oneWhistle)
                }, error = function(e) {
                    # browser()
                    cat('Error in function ', names(whistleFuns)[f],
                        ': ', as.character(e$call[1]), '-', e$message,
                        '. UID: ', oneWhistle$UID, sep='')
                }
                )
            })
        )
    }
    bind_cols(allWhistles)
}

# ceps. THIS IS ALL TERRIBLE BULLSHIT WHY ARE THEY THE SAME FUCKING FIX IT
doCepstrumCalcs <- function(cepstrumData, cepstrumFuns) {
    ### RIGHT NOW FFT NOT NEEDED - CANCELS OUT WITH INVERSE AT SAME LENGTH
    ### THIS ISNT GENERALLY TRUE, BUT IS FOR PG IMPLEMENTATION
    ### DOING IT JUST BECAUSE I CAN

    tempData <- cepstrumData[[1]]
    if(tempData$sliceData[[1]]$sliceNumber == 0) {
        tempData <- cepstrumData[[2]]
    }
    fftHop <- (tempData$startSample + 1)/tempData$sliceData[[1]]$sliceNumber
    fftLen <- tempData$sampleDuration - (tempData$nSlices - 1) * fftHop

    allCeps <- vector('list', length=length(cepstrumFuns))
    for(f in seq_along(cepstrumFuns)) {
        allCeps[[f]] <- bind_rows(
            lapply(cepstrumData, function(oneCeps) {
                oneCeps$quefrency <- oneCeps$contour
                oneCeps$time <- sapply(oneCeps$sliceData,
                                          function(x) x$sliceNumber) * fftHop / oneCeps$sr
                tryCatch({
                    cepstrumFuns[[f]](oneCeps)
                }, error = function(e) {
                    # browser()
                    cat('Error in function ', names(cepstrumFuns)[f],
                        ': ', as.character(e$call[1]), '-', e$message,
                        '. UID: ', oneCeps$UID, sep='')
                }
                )
            })
        )
    }
    bind_cols(allCeps)
}
