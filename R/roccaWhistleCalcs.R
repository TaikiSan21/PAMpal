#' @title Calculate a Set of Measurements for Whistles
#'
#' @description Calculate a set of measurements from a whistle contour. All
#'   calculations following ROCCA method from Julie and Michael Oswald, as
#'   implemented in Pamguard.
#'
#' @param data a list that must have \code{freq} the whistle contour stored as a
#'   vector of FFT bin frequencies in hertz, and \code{time} the time in seconds
#'   at each bin.
#'
#' @return A list with 50 calculated rocca parameters, each item in the list will
#'   only have 1 entry so that this can easily be converted to a data frame.
#'
#' @author Taiki Sakai \email{taiki.sakai@@noaa.gov}
#'
#' @importFrom stats sd
#' @export
#'
roccaWhistleCalcs <- function(data) {
    # ALL SLOPE CALCS ARE VERY SLIGHTLY OFF???? NOT SURE HOW TO TEST
    # I think slope difference is from getTime() getting more accuracte /
    # granular times than just step size / SR.
    # Sweep frequency calcs are off by the sweeps set to default  in te
    # setSweep() function, is the Java doing what it should?? those shouldnt
    # be counted as sweeps of default type???
    neededVals <- c('freq', 'time')
    missingVals <- neededVals[!(neededVals %in% names(data))]
    if(length(missingVals) > 0) {
        stop('Values for', paste(missingVals, collapse=', '), 'are missing.',
             'These are required for Rocca Whistle Calculations, please fix.')
    }
    contour <- data.frame(freq = data$freq, time = data$time)
    nSlices <- nrow(contour)

    result <- list(freqBeg = contour$freq[1],
                         freqEnd = contour$freq[nSlices],
                         freqMean = mean(contour$freq),
                         freqStdDev = sd(contour$freq),
                         duration = (contour$time[nSlices] - contour$time[1]))
    contour$sweep <- 0L # 0 1 2 down flat up
    contour$step <- 0L # 0 1 2 flat up down
    contour$slope <- 0
    contour$timediff <- 0
    contour$timediff[2:nSlices] <- contour$time[2:nSlices] - contour$time[1:(nSlices-1)]

    slopeBad <- contour$timediff == 0
    contour$slope[slopeBad] <- 0
    contour$slope[!slopeBad] <- (contour$freq[2:nSlices] - contour$freq[1:(nSlices-1)]) /
        (contour$timediff[2:nSlices])
    ##### FIXING HERE #####
    contour$sweep <- setSweep(contour$freq)
    inflTimeArray <- rep(NA, nrow(contour))
    for(i in 2:nSlices) {
        contour$step[i] <- setStep(contour, i, stepSens = 11)
        # Inflection / Direction calcs Java 608-676
        if(i %in% 2:3) {
            # init at 3, 1st is not real sweep so compare to prev starts at i==3
            direction <- contour$sweep[2]
        } else if((direction == 0L & contour$sweep[i] == 2L) ||
                  (direction == 2L & contour$sweep[i] == 0L)) {
            inflTimeArray[i] <- contour$time[i]
            direction <- contour$sweep[i]
        } else if(direction == 1L) {
            direction <- contour$sweep[i]
        }
    }

    # Not first, not if time diff==0. First is 0. Java 749
    useSlope <- contour$slope[contour$timediff != 0]
    result$freqSlopeMean <- mean(useSlope)
    result$freqAbsSlopeMean <- mean(abs(useSlope))
    # positive slopes
    useSlope <- contour$slope[contour$slope > 0]
    result$freqPosSlopeMean <- ifelse(length(useSlope)==0, 0,
                                      mean(useSlope))
    # negative slopes
    useSlope <- contour$slope[contour$slope < 0]
    if(length(useSlope)==0) {
        result$freqNegSlopeMean <- 0
        result$freqSlopeRatio <- 0
    } else {
        result$freqNegSlopeMean <- mean(useSlope)
        result$freqSlopeRatio <- result$freqPosSlopeMean / result$freqNegSlopeMean
    }
    # This should be fine, check Java 494
    result$freqStepUp <- sum(contour$step==1L)
    result$freqStepDown <- sum(contour$step==2L)

    result <- c(result, countSweeps(contour$sweep))
    result$numInflections <- sum(!is.na(inflTimeArray))

    # Java 694
    freqCofmSum <- 0
    if(nSlices >= 7) {
        for(i in seq(7, nSlices, 3)) { # java starts at ix=6 for some reason
            freqCofmSum <- freqCofmSum + abs(contour$freq[i]-contour$freq[i-3])
        }
    }
    result$freqCofm <- freqCofmSum / 10000

    freqQuarters <- contour$freq[round(c(1, 2, 3)*nSlices/4, 0)]
    result$freqQuarter1 <- freqQuarters[1]
    result$freqQuarter2 <- freqQuarters[2]
    result$freqQuarter3 <- freqQuarters[3]

    freqSort <- sort(contour$freq)
    result$freqSpread <- diff(quantile(freqSort, c(.25, .75), names=FALSE))
    result$freqMin <- freqSort[1]
    result$freqMax <- freqSort[nSlices]
    result$freqRange <- result$freqMax - result$freqMin
    result$freqMedian <- median(freqSort)
    result$freqCenter <- result$freqMin + result$freqRange/2
    result$freqRelBw <- result$freqRange / result$freqCenter
    result$freqMaxMinRatio <- result$freqMax / result$freqMin
    result$freqBegEndRatio <- result$freqBeg / result$freqEnd

    # Java 739
    result$freqNumSteps <- result$freqStepUp + result$freqStepDown
    result$stepDur <- result$freqNumSteps / result$duration

    # Java 770
    slopeBegAve <- mean(contour$slope[2:4]) # not first cuz 0
    if(slopeBegAve > 0) {
        result$freqBegSweep <- 2
        result$freqBegUp <- TRUE
        result$freqBegDwn <- FALSE
    } else if(slopeBegAve < 0) {
        result$freqBegSweep <- 0
        result$freqBegUp <- FALSE
        result$freqBegDwn <- TRUE
    } else {
        result$freqBegSweep <- 1
        result$freqBegUp <- FALSE
        result$freqBegDwn <- FALSE
    }
    slopeEndAve <- mean(contour$slope[(nSlices-2):nSlices])
    if(slopeEndAve > 0) {
        result$freqEndSweep <- 2
        result$freqEndUp <- TRUE
        result$freqEndDwn <- FALSE
    } else if(slopeEndAve < 0) {
        result$freqEndSweep <- 0
        result$freqEndUp <- FALSE
        result$freqEndDwn <- TRUE
    } else {
        result$freqEndSweep <- 1
        result$freqEndUp <- FALSE
        result$freqEndDwn <- FALSE
    }

    # java 810
    useSweep <- contour$sweep[2:(nSlices-1)] # not first or last, they default to 0
    if(length(useSweep) > 0) {
        result$freqSweepUpPercent <- sum(useSweep == 2L) / length(useSweep) * 100
        result$freqSweepDwnPercent <- sum(useSweep == 0L) / length(useSweep) * 100
        result$freqSweepFlatPercent <- sum(useSweep == 1L) / length(useSweep) * 100
    } else { # just in case we dont have any, set to 0
        result$freqSweepUpPercent <- 0
        result$freqSweepDwnPercent <- 0
        result$freqSweepFlatPercent <- 0
    }
    # Java 857
    if(result$numInflections > 1) {
        inflOnly <- inflTimeArray[!is.na(inflTimeArray)]
        lastInfl <- inflOnly[result$numInflections]
        smallArray <- (lastInfl - inflOnly[-result$numInflections])
        smallArray <- sort(smallArray)
        result$inflMaxDelta <- smallArray[length(smallArray)]
        result$inflMinDelta <- smallArray[1]
        result$inflMaxMinDelta <- result$inflMaxDelta / result$inflMinDelta
        result$inflMeanDelta <- mean(smallArray)
        result$inflStdDevDelta <- ifelse(length(smallArray) > 1, sd(smallArray), 0) # 0 not NA for 1 ting
        result$inflMedianDelta <- median(smallArray)
        result$inflDur <- result$numInflections / result$duration
    } else {
        result$inflMaxDelta <- 0
        result$inflMinDelta <- 0
        result$inflMaxMinDelta <- 0
        result$inflMeanDelta <- 0
        result$inflStdDevDelta <- 0
        result$inflMedianDelta <- 0
        result$inflDur <- 0
    }
    result
}

# I dont think Java comment is correct on steps covering multiple time steps 494
# This always resets to flat
setStep <- function(contour, ix, stepSens = 11) {
    # 0 1 2 flat up down
    if(ix == 1) {
        return(0L)
    }
    if(contour$step[ix-1] == 0 &&
       contour$freq[ix] >= contour$freq[ix-1]*(1+stepSens/100)) {
        return(1L)
    }
    if(contour$step[ix-1] == 0 &&
       contour$freq[ix] <= contour$freq[ix-1]*(1-stepSens/100)) {
        return(2L)
    }
    0L
}

setSweep <- function(freq) {
    nFreq <- length(freq)
    freqBefore <- freq[1:(nFreq-2)]
    freqNow <- freq[2:(nFreq-1)]
    freqAfter <- freq[3:nFreq]
    freqSweep <- rep(NA_integer_, nFreq-2)
    # 0 1 2 down flat up
    freqSweep[(freqBefore >= freqNow) & (freqNow >= freqAfter)] <- 0L
    freqSweep[(freqBefore <= freqNow) & (freqNow <= freqAfter)] <- 2L
    freqSweep[(freqBefore == freqNow) & (freqNow == freqAfter)] <- 1L
    # We've only set the 'middle' sweeps, ends don't have valid values so default to 0.
    freqSweep <- c(1L, freqSweep, 1L)
    for(i in seq_along(freqSweep)) {
        if(is.na(freqSweep[i])) {
            freqSweep[i] <- freqSweep[i-1]
        }
    }
    freqSweep
}

countSweeps <- function(sweep) {
    # 0 1 2 dwn flat up
    nSlices <- length(sweep)
    swpBefore <- sweep[2:(nSlices-1)]
    swpAfter <- sweep[3:nSlices]
    result <- list(numSweepsDwnFlat = sum(swpBefore == 0L & swpAfter == 1L))
    result$numSweepsDwnUp <- sum(swpBefore == 0L & swpAfter == 2L)
    result$numSweepsFlatDwn <- sum(swpBefore == 1L & swpAfter == 0L)
    result$numSweepsFlatUp <- sum(swpBefore == 1L & swpAfter == 2L)
    result$numSweepsUpDwn <- sum(swpBefore == 2L & swpAfter == 0L)
    result$numSweepsUpFlat <- sum(swpBefore == 2L & swpAfter == 1L)

    result
}
