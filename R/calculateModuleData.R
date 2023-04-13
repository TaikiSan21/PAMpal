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
#' @param settings a list of settings from a Pamguard XML file
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
#'
calculateModuleData <- function(binData, binFuns=list('ClickDetector'=list(standardClickCalcs)), settings=NULL) {
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
       grepl('cepstrum|burst pulse|burstpulse', detName, ignore.case = TRUE)) {
        moduleType <- 'Cepstrum'
    }
    detName <- gsub(' ', '_', detName)
    detSettings <- settings$detectors[[detName]]
    # if(!(moduleType %in% names(binFuns)) ||
    #    length(binFuns[[moduleType]])==0) {
    #     # warning("I don't have functions for Module Type ", moduleType)
    #     # If nothing, just UID and detectorName? Fine for now
    #     result <- data.frame(UID = sapply(binData$data, function(x) x$UID),
    #                          UTC = sapply(binData$data, function(x) x$date))
    #     result$UID <- as.character(result$UID)
    #     result$UTC <- convertPgDate(result$UTC)
    #     result$detectorName <- rep(detName, nrow(result))
    #     result$BinaryFile <- rep(binData$fileInfo$fileName, nrow(result))
    #     return(result)
    # }
    # Adding this to binFuns for PG version so we always at least get UID and time
    # ClickDetector gets 1 row for each channel, has 'nChan' in it
    # getBasic <- function(type) {
    #     switch(type,
    #            'ClickDetector' = function(x) {
    #                nRep <- ncol(x$wave)
    #                chan <- 1:ncol(x$wave)
    #                if('channelMap' %in% names(x)) {
    #                    mapChans <- cmapToChan(x$channelMap)
    #                    if(length(mapChans) == nRep) {
    #                        chan <- mapChans
    #                    }
    #                }
    #
    #                result <- data.frame(UID = rep(as.character(x$UID), nRep),
    #                                     UTC = rep(x$date, nRep),
    #                                     stringsAsFactors = FALSE)
    #                result$Channel <- as.character(chan)
    #                angs <- x$angles
    #                if(length(angs) == 2) {
    #                    pamWarning('More than one angle for click UID' , result$UID[1],
    #                            ', only using first.')
    #                }
    #                result$angle <- angs[1]
    #                result$angleError <- x$angleErrors[1]
    #                result
    #            }, function(x) {
    #                nRep <- if('nChan' %in% names(x)) {
    #                    x$nChan
    #                } else {
    #                    1
    #                }
    #                data.frame(UID = rep(as.character(x$UID), nRep),
    #                           UTC = rep(x$date, nRep),
    #                           stringsAsFactors = FALSE)
    #            })
    #
    # }
    # if you add new modules need to add to addBinaries - weve filtered them out there
    result <- switch(
        moduleType,
        'ClickDetector' = {
            # allClicks <- doClickCalcs(binData$data, c(getBasic(moduleType), binFuns[['ClickDetector']]), detSettings)
            allClicks <- doCalcs(binData$data, binFuns, detSettings, module='ClickDetector')
            # We want each 'type' of click to be separate 'detector'. This is PG only.
            allNames <- bind_rows(
                lapply(binData$data[unique(as.character(allClicks$UID))], function(x) {
                    data.frame(UID=as.character(x$UID),
                               detectorName=unique(c(x$type, unlist(x$annotations$classification))),
                               stringsAsFactors = FALSE)
                })) %>%
                mutate(detectorName = paste(detName, .data$detectorName, sep='_'))
            allClicks <- left_join(allClicks, allNames, join_by('UID'=='UID'), multiple='all')
            allClicks$callType <- 'click'
            allClicks
        },
        'WhistlesMoans' = {
            # allWhistles <- doWhistleCalcs(binData$data, c(getBasic(moduleType), binFuns[['WhistlesMoans']]))
            allWhistles <- doCalcs(binData$data, binFuns, detSettings, module='WhistlesMoans')
            allWhistles$detectorName <- detName
            allWhistles$callType <- 'whistle'
            allWhistles
        },
        'Cepstrum' = {
            # allCepstrum <- doCepstrumCalcs(binData$data, c(getBasic(moduleType), binFuns[['Cepstrum']]), detSettings)
            allCepstrum <- doCalcs(binData$data, binFuns, detSettings, module='Cepstrum')
            allCepstrum$detectorName <- detName
            allCepstrum$callType <- 'cepstrum'
            allCepstrum
        },
        'GPLDetector' = {
            # allGPL <- doGPLCalcs(binData$data, c(getBasic(moduleType), binFuns[['GPLDetector']]))
            allGPL <- doCalcs(binData$data, binFuns, module='GPLDetector')
            if(is.null(allGPL) ||
               nrow(allGPL) == 0) {
                return(NULL)
            }
            allGPL$detectorName <- detName
            allGPL$callType <- 'gpl'
            allGPL
        },
        pamWarning("I don't know how to deal with Module Type ", moduleType)
    )
    result$BinaryFile <- binData$fileInfo$fileName
    result$UTC <- convertPgDate(result$UTC)
    result
}

doCalcs <- function(data, funs, detSettings=NULL, module, retry=TRUE) {
    getBasic <- function(type) {
        switch(type,
               'ClickDetector' = function(x) {
                   nRep <- ncol(x$wave)
                   chan <- 1:ncol(x$wave)
                   if('channelMap' %in% names(x)) {
                       mapChans <- cmapToChan(x$channelMap)
                       if(length(mapChans) == nRep) {
                           chan <- mapChans
                       }
                   }

                   result <- data.frame(UID = rep(as.character(x$UID), nRep),
                                        UTC = rep(x$date, nRep),
                                        stringsAsFactors = FALSE)
                   result$Channel <- as.character(chan)
                   angs <- x$angles
                   if(length(angs) == 2) {
                       pamWarning('More than one angle for click UID' , result$UID[1],
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
    funs <- c(getBasic(module), funs[[module]])
    names(funs)[1] <- 'getBasic'
    switch(module,
           'ClickDetector' = {
               if(isTRUE(detSettings$decimated)) {
                   for(i in seq_along(data)) {
                       data[[i]]$sr <- detSettings$sr
                   }
                   for(f in seq_along(funs)) {
                       if('sr_hz' %in% names(formals(funs[[f]]))) {
                           formals(funs[[f]])[['sr_hz']] <- 'auto'
                       }
                   }
               }
           },
           'WhistlesMoans' = {
               if(is.null(detSettings)) {
                   data <- contourToFreq(data)
               } else {
                   for (i in seq_along(data)) {
                       # this not perfect for changing SR so we should just not do it for whistles
                       sr <- detSettings$sr
                       fftLen <- detSettings$length
                       fftHop <- detSettings$hop
                       data[[i]]$freq <- data[[i]]$contour * sr/fftLen
                       data[[i]]$allFreq <- do.call(cbind, lapply(data[[i]]$sliceData,
                                                                         function(x) x$peakData)) * sr/fftLen
                       data[[i]]$time <- sapply(data[[i]]$sliceData, function(x) x$sliceNumber) *
                           fftHop/sr
                   }
               }
           },
           'Cepstrum' = {
               if(is.null(detSettings)) {
                   tempData <- data[[1]]
                   if(tempData$sliceData[[1]]$sliceNumber == 0) {
                       tempData <- data[[2]]
                   }
                   fftHop <- (tempData$startSample + 1)/tempData$sliceData[[1]]$sliceNumber
                   fftLen <- tempData$sampleDuration - (tempData$nSlices - 1) * fftHop
               } else {
                   fftHop <- detSettings$hop
                   fftLen <- detSettings$length
                   if(isTRUE(detSettings$decimated)) {
                       for(i in seq_along(data)) {
                           data[[i]]$sr <- detSettings$sr
                       }
                   }
               }
               for(i in seq_along(data)) {
                   data[[i]]$quefrency <- data[[i]]$contour / data[[i]]$sr
                   data[[i]]$time <- sapply(data[[i]]$sliceData, function(s) s$sliceNumber) * fftHop / data[[i]]$sr
               }
           },
           'GPLDetector' = {
               if(is.null(detSettings)) {
                   freqRes <- data[[1]]$freqRes
                   timeRes <- data[[1]]$timeRes
               } else {
                   # this not perfect for changing SR so we should just not do it for whistles
                   sr <- detSettings$sr
                   fftLen <- detSettings$length
                   fftHop <- detSettings$hop
                   freqRes <- sr / fftLen
                   timeRes <- fftHop / sr
                   # whistleData[[i]]$freq <- whistleData[[i]]$contour * sr/fftLen
                   # whistleData[[i]]$allFreq <- do.call(cbind, lapply(whistleData[[i]]$sliceData,
                   #                                                   function(x) x$peakData)) * sr/fftLen
                   # whistleData[[i]]$time <- sapply(whistleData[[i]]$sliceData, function(x) x$sliceNumber) *
                   #     fftHop/sr

               }

               for(i in seq_along(data)) {
                   # if(gplData[[i]]$area <= 2) {
                   #     next
                   # }
                   allTimes <- data[[i]]$points[1, ]
                   times <- unique(allTimes)
                   freq <- rep(NA, length(times))
                   thisEnergy <- data[[i]]$energy
                   thisFreq <- data[[i]]$points[2, ] * freqRes
                   for(j in seq_along(times)) {
                       whichThis <- allTimes == times[j]
                       maxEn <- max(thisEnergy[whichThis])
                       whichMax <- whichThis & (thisEnergy == maxEn)
                       freq[j] <- thisFreq[whichMax[1]]
                   }
                   data[[i]]$freq <- freq
                   data[[i]]$time <- times * timeRes
                   if(length(times) == 1) {
                       data[[i]]$timeRes <- timeRes
                   }
               }
           }
    )
    allDets <- vector('list', length=length(funs))
    for(f in seq_along(funs)) {
        allDets[[f]] <- bind_rows(
            lapply(data, function(oneDet) {
                tryCatch({
                    funs[[f]](oneDet)
                }, error = function(e) {
                    if(isFALSE(retry)) {
                        message('Error in function ', names(funs)[f],
                                ':', as.character(e$call[1]), '-', e$message,
                                '. UID: ', oneDet$UID, sep='')
                    }
                })
            })
        )
    }
    nDets <- sapply(allDets, nrow)
    if(isTRUE(retry) &&
       !all(nDets == nDets[1])) {
        return(doCalcs(data, funs, detSettings, module, retry=FALSE))
    }
    bind_cols(allDets)
}



# In general just make sure data has $wave and $sr for clicks?
# doClickCalcs <- function(clickData, clickFuns, detSettings=NULL) {
#     allClicks <- vector('list', length = length(clickFuns))
#     # This just returns a df we will bind to db by UID
#     if(isTRUE(detSettings$decimated)) {
#         for(i in seq_along(clickData)) {
#             clickData[[i]]$sr <- detSettings$sr
#         }
#         for(f in seq_along(clickFuns)) {
#             if('sr_hz' %in% names(formals(clickFuns[[f]]))) {
#                 formals(clickFuns[[f]])[['sr_hz']] <- 'auto'
#             }
#         }
#     }
#     for(f in seq_along(clickFuns)) {
#         # Apply this function to each datapoint in binary
#         allClicks[[f]] <- bind_rows(
#             lapply(clickData, function(oneClick) {
#                 tryCatch({
#                     clickFuns[[f]](oneClick)
#                 }, error = function(e) {
#                     # browser()
#                     message('Error in function ', names(clickFuns)[f],
#                         ': ', as.character(e$call[1]), '-', e$message,
#                         '. UID: ', oneClick$UID, sep='')
#                 }
#                 )
#             })
#         )
#     }
#     bind_cols(allClicks)
# }
#
#
# # WHISTLES DAKINE
# # Hop = (StartSample + 1) / (1st Slice Num)
# # FFT Length = Samp Dur - (N slices - 1) * Hop(above)
# # FFT Length = SampleRate * max(peakData) / MaxFreq
# # Mult bin by SR / Len
#
# # PROBLEM right now this is very PG specific. Which is maybe okay, because
# # moduleCalcs is PG-based? Possbly rename these doPgWhistleCalcs or w/e
# doWhistleCalcs <- function(whistleData, whistleFuns, detSettings=NULL) {
#     # Does this die if 1st slice is #1 (or 0, w/e index is)
#     # probably, so use next whistle
#     if(is.null(detSettings)) {
#         whistleData <- contourToFreq(whistleData)
#     } else {
#         for (i in seq_along(whistleData)) {
#             # this not perfect for changing SR so we should just not do it for whistles
#             sr <- detSettings$sr
#             fftLen <- detSettings$length
#             fftHop <- detSettings$hop
#             whistleData[[i]]$freq <- whistleData[[i]]$contour * sr/fftLen
#             whistleData[[i]]$allFreq <- do.call(cbind, lapply(whistleData[[i]]$sliceData,
#                                                        function(x) x$peakData)) * sr/fftLen
#             whistleData[[i]]$time <- sapply(whistleData[[i]]$sliceData, function(x) x$sliceNumber) *
#                 fftHop/sr
#         }
#     }
#     # tempData <- whistleData[[1]]
#     # if(tempData$sliceData[[1]]$sliceNumber == 0) {
#     #     tempData <- whistleData[[2]]
#     # }
#     # fftHop <- (tempData$startSample + 1)/tempData$sliceData[[1]]$sliceNumber
#     # fftLen <- tempData$sampleDuration -
#     #     (tempData$sliceData[[tempData$nSlices]]$sliceNumber - tempData$sliceData[[1]]$sliceNumber) * fftHop
#
#     allWhistles <- vector('list', length=length(whistleFuns))
#     for(f in seq_along(whistleFuns)) {
#         allWhistles[[f]] <- bind_rows(
#             lapply(whistleData, function(oneWhistle) {
#                 tryCatch({
#                     whistleFuns[[f]](oneWhistle)
#                 }, error = function(e) {
#                     # browser()
#                     message('Error in function ', names(whistleFuns)[f],
#                         ': ', as.character(e$call[1]), '-', e$message,
#                         '. UID: ', oneWhistle$UID, sep='')
#                 }
#                 )
#             })
#         )
#     }
#     bind_cols(allWhistles)
# }
#
# # ceps. THIS IS ALL TERRIBLE BULLSHIT WHY ARE THEY THE SAME FUCKING FIX IT
# doCepstrumCalcs <- function(cepstrumData, cepstrumFuns, detSettings=NULL) {
#     ### RIGHT NOW FFT NOT NEEDED - CANCELS OUT WITH INVERSE AT SAME LENGTH
#     ### THIS ISNT GENERALLY TRUE, BUT IS FOR PG IMPLEMENTATION
#     ### DOING IT JUST BECAUSE I CAN
#     if(is.null(detSettings)) {
#         tempData <- cepstrumData[[1]]
#         if(tempData$sliceData[[1]]$sliceNumber == 0) {
#             tempData <- cepstrumData[[2]]
#         }
#         fftHop <- (tempData$startSample + 1)/tempData$sliceData[[1]]$sliceNumber
#         fftLen <- tempData$sampleDuration - (tempData$nSlices - 1) * fftHop
#     } else {
#         fftHop <- detSettings$hop
#         fftLen <- detSettings$length
#         if(isTRUE(detSettings$decimated)) {
#             for(i in seq_along(cepstrumData)) {
#                 cepstrumData[[i]]$sr <- detSettings$sr
#             }
#         }
#     }
#
#     allCeps <- vector('list', length=length(cepstrumFuns))
#     for(f in seq_along(cepstrumFuns)) {
#         allCeps[[f]] <- bind_rows(
#             lapply(cepstrumData, function(oneCeps) {
#                 oneCeps$quefrency <- oneCeps$contour
#                 oneCeps$time <- sapply(oneCeps$sliceData,
#                                           function(x) x$sliceNumber) * fftHop / oneCeps$sr
#                 tryCatch({
#                     cepstrumFuns[[f]](oneCeps)
#                 }, error = function(e) {
#                     # browser()
#                     message('Error in function ', names(cepstrumFuns)[f],
#                         ': ', as.character(e$call[1]), '-', e$message,
#                         '. UID: ', oneCeps$UID, sep='')
#                 }
#                 )
#             })
#         )
#     }
#     bind_cols(allCeps)
# }
#
# # GPL
# doGPLCalcs <- function(gplData, gplFuns, detSettings=NULL) {
#     # Does this die if 1st slice is #1 (or 0, w/e index is)
#     # probably, so use next whistle
#     if(is.null(detSettings)) {
#         freqRes <- gplData[[1]]$freqRes
#         timeRes <- gplData[[1]]$timeRes
#     } else {
#         # this not perfect for changing SR so we should just not do it for whistles
#         sr <- detSettings$sr
#         fftLen <- detSettings$length
#         fftHop <- detSettings$hop
#         freqRes <- sr / fftLen
#         timeRes <- fftHop / sr
#         # whistleData[[i]]$freq <- whistleData[[i]]$contour * sr/fftLen
#         # whistleData[[i]]$allFreq <- do.call(cbind, lapply(whistleData[[i]]$sliceData,
#         #                                                   function(x) x$peakData)) * sr/fftLen
#         # whistleData[[i]]$time <- sapply(whistleData[[i]]$sliceData, function(x) x$sliceNumber) *
#         #     fftHop/sr
#
#     }
#
#     for(i in seq_along(gplData)) {
#         # if(gplData[[i]]$area <= 2) {
#         #     next
#         # }
#         allTimes <- gplData[[i]]$points[1, ]
#         times <- unique(allTimes)
#         freq <- rep(NA, length(times))
#         thisEnergy <- gplData[[i]]$energy
#         thisFreq <- gplData[[i]]$points[2, ] * freqRes
#         for(j in seq_along(times)) {
#             whichThis <- allTimes == times[j]
#             maxEn <- max(thisEnergy[whichThis])
#             whichMax <- whichThis & (thisEnergy == maxEn)
#             freq[j] <- thisFreq[whichMax]
#         }
#         gplData[[i]]$freq <- freq
#         gplData[[i]]$time <- times * timeRes
#         if(length(times) == 1) {
#             gplData[[i]]$timeRes <- timeRes
#         }
#     }
#     # tempData <- whistleData[[1]]
#     # if(tempData$sliceData[[1]]$sliceNumber == 0) {
#     #     tempData <- whistleData[[2]]
#     # }
#     # fftHop <- (tempData$startSample + 1)/tempData$sliceData[[1]]$sliceNumber
#     # fftLen <- tempData$sampleDuration -
#     #     (tempData$sliceData[[tempData$nSlices]]$sliceNumber - tempData$sliceData[[1]]$sliceNumber) * fftHop
#     allGpl <- vector('list', length=length(gplFuns))
#     for(f in seq_along(gplFuns)) {
#         allGpl[[f]] <- bind_rows(
#             lapply(gplData, function(oneGpl) {
#                 # if(length(oneGpl$freq) <= 2) {
#                 #     return(NULL)
#                 # }
#                 tryCatch({
#                     gplFuns[[f]](oneGpl)
#                 }, error = function(e) {
#                     # browser()
#                     message('Error in function ', names(gplFuns)[f],
#                             ': ', as.character(e$call[1]), '-', e$message,
#                             '. UID: ', oneGpl$UID, sep='')
#                 }
#                 )
#             })
#         )
#     }
#     bind_cols(allGpl)
# }
