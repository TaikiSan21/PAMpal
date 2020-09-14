#' @title Calculate a Set of Measurements for Clicks EG Style
#'
#' @description Calculate clicks per EGs methods for comparison to old
#'   clustering study. DO NOT USE THESE, ONLY CREATED FOR TESTING AND
#'   NOT EXPECTED TO BE USEFUL FOR GENERAL PUBLIC
#'
#' @param data a list that must have 'wave' containing the wave form as a
#'   matrix with a separate column for each channel, and 'sr' the
#'   sample rate of the data. Data can also be a \code{Wave} class
#'   object, like one created by \code{\link[tuneR]{Wave}}.
#' @param sr_hz either \code{'auto'} (default) or the numeric value of the sample
#'   rate in hertz. If \code{'auto'}, the sample rate will be read from the
#'   'sr' of \code{data}
#' @param calibration a calibration function to apply to the spectrum, must be
#'   a gam. If NULL no calibration will be applied (not recommended).
#' @param filterfrom_khz frequency in khz of highpass filter to apply, or the lower
#'   bound of a bandpass filter if \code{filterto_khz} is not \code{NULL}
#' @param filterto_khz if a bandpass filter is desired, set this as the upper bound.
#'   If only a highpass filter is desired, leave as the default \code{NULL} value
#'
#' @return A data frame with one row for each channel of click waveform.
#'   Calculates approximate noise level and click duration from the
#'   TK energy (Soldevilla JASA17), up to 3 highest peak frequencies and
#'   the 'troughs' between them (see \code{\link[PAMmisc]{peakTrough}}), and the 3
#'   and 10dB bandwidth levels and 'Q' value (see \code{\link[seewave]{Q}}).
#'
#' @author Taiki Sakai \email{taiki.sakai@@noaa.gov}
#'
#' @importFrom seewave bwfilter TKEO spec specprop fpeaks
#' @importFrom dplyr bind_rows
#' @importFrom stats median quantile loess
#' @importFrom tuneR WaveMC
#' @export
#'
egClicks <- function(data, sr_hz='auto', calibration=NULL, filterfrom_khz=90, filterto_khz=NULL) {

    if(inherits(data, 'Wave')) {
        data <- WaveMC(data)
        data <- list(wave = data@.Data, sr = data@samp.rate)
    }
    if(!is.matrix(data$wave)) {
        data$wave <- matrix(data$wave, ncol=1)
    }
    FFTsize <- 512
    result <- list()
    for(chan in 1:ncol(data$wave)) {
        if(sr_hz == 'auto') {
            sr <- data$sr
        } else {
            sr <- sr_hz
        }
        # i should be a binary, j each click
        DSc=NULL
        DSc <- data.frame(Channel = chan)
        # it appears FFF is a waveform already HPfilter at 90khz and only a 512 clip. around 103 Step3
        # filtWav90HP=foo[[i]][[j]] # make this just wave
        thisWav <- data$wave[, chan]
        filtWav90HP <- seewave::bwfilter(thisWav,f=sr,n=4,from=filterfrom_khz, bandpass=TRUE, output = "sample")

        # pk=min(which(grepl(max(filtWav90HP),filtWav90HP)))
        pk <- which.max(filtWav90HP)
        minn=pk-(FFTsize/2)
        #Make sure that your minimum sample isn't less than 0
        min=ifelse(minn<1,1,minn)
        maxx=pk+(FFTsize/2)-1
        #Make sure that your minimum sample isn't greater than 5000
        # these were 8k length clicks
        max=ifelse(maxx>5000,5000, maxx)
        # if(max > length(thisWav))
        #This should equal 512 samples
        filtRange <- checkNotOut(min, FFTsize, length(filtWav90HP))

        # filtWav90HP=filtWav90HP[min:max]
        filtWav90HP=filtWav90HP[filtRange]

        Nmin=ifelse(min <=800, 2400,10)
        Nmax=Nmin+FFTsize-1

        noiseRange <- checkNotOut(Nmin, FFTsize, length(thisWav))
        # noise <- thisWav[Nmin:Nmax]
        noise <- thisWav[noiseRange]

        # noise is a 512 clip from start or end of wav depending where click is, 124 step3
        # noise=eventnoise[[i]][[j]] # make this just wave

        #Create spectrum.
        thisSpec=seewave::spec(filtWav90HP,f=sr ,wl=FFTsize, norm = FALSE, correction = "amplitude", plot=F)
        #Convert amplitude to relative dB.
        reldB=(20*log10(thisSpec[,2]))
        #Conver -Inf to NA
        reldB[!is.finite(reldB)] <- NA

        #Calibration Curve
        newClick=data.frame(Freq=(thisSpec[,1]*1000),Sensitivity = reldB)

        #Apply GAMs from both the HP and the recorder.
        if(!is.null(calibration)) {
            #DO CA
            if(is.function(calibration)) {
                calFun <- calibration
            } else if(is.character(calibration)) {
                calFun <- findCalibration(calibration)
            }
        } else {
            calFun <- function(x) {
                0
            }
        }

        TKfunc.left= seewave::TKEO(filtWav90HP,f=sr,M=1,plot = F)

        # convert energy to dB scale (note multiplier is 10 for energy, 20 for sound pressure)
        TKenergy=TKfunc.left[,2]

        TKenergyDB= 10*log10(TKenergy-min(TKenergy,na.rm=TRUE))
        # normalize to 0 dB max
        TKenergyDB= TKenergyDB - max(TKenergyDB,na.rm=TRUE)
        TKenergyDB[!is.finite(TKenergyDB)] <- NA

        #Apply combined curve
        CaliTKenergyDB= TKenergyDB + calFun(seq(from=0, by = TKfunc.left[2,1] - TKfunc.left[1,1], length.out = nrow(TKfunc.left)) * 1e3)

        #Find noise median
        noiselevel=median(CaliTKenergyDB,na.rm=TRUE)
        #Remove clicks with a noise median greater than -15 dB.
        # if (noiselevel > -10) next

        #Duration

        #Find noise threshold by multiplying all energy above the 40% threshold by 100.
        noisethreshold=quantile(TKfunc.left[,2],probs = .40, na.rm = TRUE)*100
        #Subset the TK function for energy above the 40% threshold.

        dur=subset(TKfunc.left,TKfunc.left[,2]>= noisethreshold)
        #Subtract the max time value from the minimum time value for the click duration.
        if(length(dur) == 0) {
            duration <- 0
        } else {
            duration=1000000*(max(dur[,1])-min(dur[,1]))
        }

        #Repeat for noise sample
        nar=seewave::spec(noise,f=sr,wl=FFTsize, norm = FALSE, correction = "amplitude", plot=F)
        #Convert amplitude to relative dB.
        nreldB=(20*log10(nar[,2]))
        #Conver -Inf to NA
        nreldB[!is.finite(nreldB)] <- NA

        #Calibration Curve
        newNoise=data.frame(Freq=(nar[,1]*1000),Sensitivity = nreldB)

        #Applying the calibration curve to the click and noise sample
        clicksens=reldB + calFun(seq(from=0, by = thisSpec[2,1] - thisSpec[1,1], length.out = nrow(thisSpec)) * 1e3)
        noisesens=nreldB + calFun(seq(from=0, by = nar[2,1] - nar[1,1], length.out = nrow(nar)) * 1e3)
        #Smooth the jagged noise sensitivity for each known frequnecy measure.
        if(any(is.na(noisesens))) {
            sensdiff <- clicksens
        } else {
            x=loess(noisesens~newNoise[,1], family="gaussian")
            #Calculate the difference between the click sensitivity and the smoothed noise sample.
            sensdiff=(clicksens-x$fitted)
            #Remove extremes in the lower frequencies
            sensdiff[1:50]=0
        }

        Adj4zero=sensdiff-max(sensdiff)

        #Calibrated spectrum for futher feature measures.
        TotCalibr= cbind(newClick$Freq/1000,Adj4zero)

        #Determine if there is a shoulder present and calculate it's frequency.
        shod=shoulder(TotCalibr,f=sr,threshold = -15,amp.slope = 0.1, freq.dist = 5, plot=F)

        #Find the Q and Bandwidth data for -3,-6, and -10 dB.
        #From this function we get the peak Hz, center Hz, bandwidth, and BW max and min.
        dbBW3=Qfast(TotCalibr,f=sr,level=-3, plot = F)
        # a <- lapply(dbBW3, function(x) ifelse(is.null(x), NA, x))
        # stats3=as.data.frame(t(unlist(a)))
        names(dbBW3)=c("Q_3db", "PeakHz_3dB", "fmin_3db","fmax_3db", "BW_3db")
        dbBW3$centerHz_3db=dbBW3$fmax_3db-(dbBW3$BW_3db/2)

        dbBW6=Qfast(TotCalibr,f=sr,level=-6, plot =F)
        # b <- lapply(dbBW6, function(x) ifelse(is.null(x), NA, x))
        # stats6=as.data.frame(t(unlist(b)))
        names(dbBW6)=c("Q_6db", "PeakHz_6dB", "fmin_6db","fmax_6db", "BW_6db")
        dbBW6$centerHz_6db=dbBW6$fmax_6db-(dbBW6$BW_6db/2)

        #Bind all together.
        DSc=cbind(DSc,dbBW3,dbBW6)  #stats10,

        #find the QRMS
        stats = seewave::specprop(thisSpec,f=sr, flim=c(min(100, sr/2e3), min(144, sr/2e3)))
        #spectral standard deviation around (+/-) the centroid frequency (Kyhn et al 2013).
        qrms=DSc$centerHz_3db/(stats$sd*2)

        #Bring it all together
        DSc$duration=duration
        DSc$noiselevel=noisethreshold
        DSc$QRMS=qrms
        DSc$bwRMS=stats$sd
        shoulder.kHz=ifelse(length(shod)>0,shod[1],0)
        DSc$shoulder.kHz=ifelse(shoulder.kHz==0,0,shoulder.kHz-DSc$PeakHz_3dB)
        DSc$shoulder.amp=ifelse(length(shod)>0,shod[2],0)
        result[[chan]] <- DSc
    }
    result <- dplyr::bind_rows(result)
    result$Channel <- as.character(result$Channel)
    result
}

shoulder=function(spec,f,threshold,amp.slope,freq.dist, plot) {
    #Find peak Hz
    prim=seewave::fpeaks(spec,f=f,nmax=1, plot=F)
    #Find all peaks above your amplitude threshold.
    pks=seewave::fpeaks(spec,f=f, threshold = threshold,  amp=c(amp.slope,amp.slope), plot=F)
    #Subset those peaks to remove ones too close the actual peak.
    res=subset(pks,pks[,1] >= prim[1]+freq.dist | pks[,1] <= prim[1]-freq.dist)
    #Removes false shoulders below 100 kHz
    res1=subset(res,res[,1]>=100)
    #Find Shoulder
    mx=which.max(res1[,2])
    shod=res1[mx,]

    #Plot
    b.pks=rbind(prim,shod)


    shod

}

checkNotOut <- function(min, length=512, limit) {
    if(min <= 1) {
        return(1:(1+length-1))
    }
    if(min + length -1 >= limit) {
        return((limit-length+1):limit)
    }
    return(min:(min + length - 1))
}
