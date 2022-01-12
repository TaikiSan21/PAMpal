# Creating custom functions

PAMpal comes bundled with a "standard" set of a parameter calculations for each of
the four detectors currently supported - Click Detector, Whistle and Moan Detector,
Cepstrum Detector, and GPL Detector. These parameters are based on methods that have
been tested in the past and have proven useful in various species classification
tasks, but by no means do we think they are perfect. There is always room for new 
ideas, and one of the key focuses of PAMpal is to support other people's work and
make it easy to implement new methods. This takes a bit of work and some programming
knowledge, but otherwise fits seamlessly into the regular PAMpal processing flow
once a user's custom functions have been added to the [PAMpalSettings][pps] object
with `addFunction`. 

Of course we can't just add any function you might have to PAMpal, it must follow
some specific rules of the kind of data your function expects to process as well 
as the kind of output your function creates. Most of this is pretty straight forward,
and PAMpal will do a quick check to see if your function runs without crashing and
produces the right type of output on a test detection when you try to add it. 
Fair warning: Just because your function passes this test does not mean it will always
run smoothly on a large dataset. Acoustic datasets are massive, and within all those
detections there are likely to be some with inputs that you weren't expecting (e.g.
click clips that are all zeroes, whistle detections that are only 1-2 FFT bins long).

For PAMpal's "standard" functions, there are a few parameters for functions that need 
to be set when the function is added. This functionality extends to custom user-added 
functions as well - any parameters you include other than the defaults discussed below (usually
just `data`, maybe one more) will trigger the same kind of pop-up message asking users
to specify values for those parameters.

Here's a breakdown of the kind of data your function can expect for each detector, and
the format of the expected output. I've also included a blank template function that
you should be able to copy for each one to get you started. Just copy the template, 
rename it to something useful, and start filling in the pieces you need!

### Click Detector

The click detector is actually the most complicated of the detectors, but this is 
mostly because PAMpal's click processing module is designed to work with a 
frequency response calibration function if you supply one. This means that any
click processing functions need to know how to access and use that information.

The click module expects parameters `data` and `calibration`
in addition to any adjustable parameters you want for your function. If you don't need
to work with a calibration function you can leave this parameter (and associated code chunk)
out of your function. If you do want to be able to access a calibration function, simply leave in the
`calibration=NULL` in the function definition and the associated code below. 

`data` is a list with two items: 

  1. `data$wave` is a matrix of the values of the waveform, one column for each channel
  2. `data$sr` is the sample rate, in Hz

Your click function should do the same set of calculations for each channel, and
return a dataframe with 1 row for each channel. It is not necessary to include
a "Channel" column keeping track of which row is which channel, PAMpal handles this
for you elsewhere.

Here's the template with some comments and examples. **IMPORTANT** note the
`packageList` section at the top - PAMpal is designed to make it easy to share
these kinds of functions with others so that everyone can make use of your smart 
ideas, this `packageList` checks to make sure that other users have all the necessary
packages they need to run for your function (or installs if missing). Make sure to fill
this out completely to avoid headaches.

```{r}
sampleClickFun <- function(data, calibration=NULL) {
  # List names of all required packages here to guarantee they will be 
  # installed when others try to use your function, ex:
  # packageList <- c('seewave', 'signal')
  packageList <- c()
  for(p in packageList) {
    if(!require(p)) {
      install.packages(p)
      require(p)
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
```

### Whistle and Moan & GPL

Functions for the Whistle & Moan detector or GPL detector both work on the same
kind of contour data, so have the same requirements. These types of detections do
not have any calibration to worry about, so the only expected parameter is `data`.

`data` is a list with two items:

  1. `freq` a vector of frequency bins making up the contour, units in Hertz
  2. `time` a vector of the time at each bin, units in seconds
  
The expected output for this function is a dataframe with a single row or a 
named list where all items have length 1.

Here's the template, again be sure to fill in the `packageList`:
  
```{r}
sampleContourFun <- function(data) {
  # List names of all required packages here to guarantee they will be 
  # installed when others try to use your function, ex:
  # packageList <- c('seewave', 'signal')
  packageList <- c()
  for(p in packageList) {
    if(!require(p)) {
      install.packages(p)
      require(p)
    }
  }
  # prepping output
  result <- list()
  ### DO STUFF AFTER HERE ####
  # Store results in result:
  # maxFreq <- max(data$freq)
  # result$maxFreq <- maxFreq
  result
}
```

### Cepstrum

Functions for the Cepstrum module are similarly straight forward. This detector
works using the Whistle & Moan detector module, so the data structure is similar.
Just like the contour detectors above, the only expected parameter is `data`.

`data` is a list with three items:

  1. `quefrency` a vector of the quefrency bins making up the contour. Quefrency
  is the result of taking the "cepstrum" of a signal. For practical purposes it is
  enough to know that quefrency / sample rate = ICI.
  2. `time` a vector of the time at each bin, units in seconds
  3. `sr` the sample rate, units Hz
  
The expected output is a dataframe with 1 row or a named list where all items have
length 1.

```{r}
sampleCepstrumFun <- function(data) {
  # List names of all required packages here to guarantee they will be 
  # installed when others try to use your function, ex:
  # packageList <- c('seewave', 'signal')
  packageList <- c()
  for(p in packageList) {
    if(!require(p)) {
      install.packages(p)
      require(p)
    }
  }
  # prepping output
  result <- list()
  ### DO STUFF AFTER HERE ####
  # Store results in result:
  # maxICI <- max(data$quefrency / data$sr)
  # result$maxICI <- maxICI
  result
}
```

[pps]: PAMpalSettings.md