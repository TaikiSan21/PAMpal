#' @title Create Wav Clips of Data
#'
#' @description Creates audio clips containing sounds from events or detections
#'
#' @details \code{parseEventClipName} parses the file names created to pull out
#'   event names or file start times
#'
#' @param x \linkS4class{AcousticStudy} object containing data to make wav clips for
#' @param buffer amount before and after each event to also include in the clip, in seconds.
#'   Can either be a vector of length two specifying how much to buffer before and after
#'   (first number should be negative), or a single value if the buffer amount should be identical.
#' @param outDir directory to write clips to, defaults to current directory
#' @param mode either \code{'event'} or \code{'detection'} specifying whether to create
#'   wav clips of entire events or individual detections
#' @param channel channel(s) of clips to write
#' @param filter filter to apply to wav clips before writing, values in kHz. A value of \code{0}
#'   applies no filter. A single value applies a highpass filter at that value. A vector of two
#'   values applies a lowpass filter if the first number is \code{0}, or a bandpass filter if
#'   both are non-zero.
#' @param useSample logical flag to use startSample information in binaries instead of UTC
#'   time for start of detections. This can be slightly more accurate (~1ms) but will take
#'   longer
#' @param fixLength logical flag to fix the output clip length to a constant value. If
#'   \code{TRUE}, then output clip length is entirely determined by the buffer value, as
#'   if the detection or event had zero length. E.g. \code{buffer=c(-2,1)} will produce clips
#'   3 seconds long, starting 2 seconds before the detection/event start time.
#' @param progress logical flag to show progress bar
#' @param verbose logical flag to show summary messages
#'
#' @return A vector of file names for the wav clips that were successfully
#'   created, any that were not able to be written will be \code{NA}. Note
#'   that currently this can only write clips with up to 2 channels. File names
#'   will be formatted as
#'   [Event or Detection]_[Id]CH[ChannelNumber(s)]_[YYYYMMDD]_[HHMMSS]_[mmm].wav
#'   The last numbers are the start time of the file in UTC, accurate to milliseconds.
#'   The Id is either the event ID or the detection UID.
#'
#' @examples
#'
#' data(exStudy)
#' recs <- system.file('extdata', 'Recordings', package='PAMpal')
#' exStudy <- addRecordings(exStudy, folder=recs, log=FALSE, progress=FALSE)
#' \dontrun{
#' # not running so that no wav clips are written to disk
#' wavs <- writeEventClips(exStudy, outDir='WavFolder', mode='event')
#' }
#'
#' @author Taiki Sakai \email{taiki.sakai@@noaa.gov}
#'
#' @importFrom dplyr bind_rows arrange group_by summarise ungroup
#' @importFrom tuneR readWave writeWave MCnames bind
#'
#' @export
#'
writeEventClips <- function(x, buffer = c(0, 0.1), outDir='.', mode=c('event', 'detection'),
                            channel = 1, filter=0, useSample=FALSE, progress=TRUE, verbose=TRUE,
                            fixLength=FALSE) {
    if(!dir.exists(outDir)) dir.create(outDir)
    if(length(channel) > 2) {  #### WAV CLIP SPECIFIC
        message('R can only write wav files with 2 or less channels, channels will be split',
                ' across multiple files.')
        allFiles <- character(0)
        for(i in 1:(ceiling(length(channel)/2))) {
            ix <- (i-1)*2 +1
            if((ix+1) <= length(channel)) {
                thisChan <- channel[ix:(ix+1)]
            } else {
                thisChan <- channel[ix]
            }
            allFiles <- c(allFiles, writeEventClips(x, buffer=buffer, outDir=outDir, mode=mode, filter=filter,
                                                    channel = thisChan, useSample=useSample,progress=progress,
                                                    verbose=verbose, fixLength=fixLength))
        }
        return(allFiles)
    }
    getClipData(x, buffer=buffer, mode=mode, channel=channel, useSample=useSample,
                progress=progress, verbose=verbose, FUN=writeOneClip, outDir=outDir,
                filter=filter, fixLength=fixLength)
}

writeOneClip <- function(wav, name, time, channel, mode, outDir='.', filter) {
    fileName <- paste0(oneUpper(mode), '_', name, 'CH', paste0(channel, collapse=''))
    fileName <- paste0(fileName, '_',psxToChar(time[1]))
    fileName <- paste0(gsub('\\.wav$', '', fileName), '.wav')
    # timeRange[1] is actual start time in posix
    fileName <- file.path(outDir, fileName)
    if(length(filter) == 1) {
        filterFrom <- filter * 1e3
        filterTo <- NULL
    } else {
        filterFrom <- filter[1] *1e3
        filterTo <- filter[2] * 1e3
    }
    if(filterFrom == 0) {
        filterFrom <- NULL
    }
    if(!is.null(filterFrom) ||
       !is.null(filterTo)) {
        for(i in 1:ncol(wav@.Data)) {
            wav@.Data[, i] <- round(seewave::bwfilter(wav@.Data[, i], f=wav@samp.rate, from=filterFrom, to=filterTo)[, 1], 0)
        }
    }
    writeWave(wav, fileName, extensible = FALSE)
    fileName
}



oneUpper <- function(x) {
    paste0(toupper(substr(x, 1, 1)),
           tolower(substr(x, 2, nchar(x))))
}

psxToChar <- function(x) {
    psFloor <- as.character(as.POSIXct(floor(as.numeric(x)), origin='1970-01-01 00:00:00', tz='UTC'))
    psMilli <- round(as.numeric(x)-floor(as.numeric(x)), 3)
    psMilli <- sprintf('%.3f',psMilli)
    psMilli <- substr(psMilli, 3, 5)
    psFloor <- gsub('-| |:', '', psFloor)
    paste0(substr(psFloor, 1, 8), '_',
           substr(psFloor, 9, 14), '_',
           gsub('^0\\.', '', psMilli))
}

#' @rdname writeEventClips
#'
#' @param file file name to parse
#' @param part part of file name to return
#'
#' @export
#'
parseEventClipName <- function(file, part=c('event', 'time', 'UID', 'channel', 'UTC')) {
    if(length(file) > 1) {
        return(sapply(file, function(x) {
            parseEventClipName(x, part=part)
        }, USE.NAMES=FALSE))
    }
    file <- basename(file)
    pattern <- '(Event|Detection)_(.*)(CH[0-9]{1,2})_([0-9]{14}_[0-9]{3}|[0-9]{8}_[0-9]{6}_[0-9]{3})\\.wav$'
    switch(match.arg(part),
           'event' = {
               result <- gsub(pattern, '\\2', file)
               result <- strsplit(result, '\\.')[[1]]
               result <- paste0(result[1], '.', result[2])
               result
           },
           'UID' = {
               if(gsub(pattern, '\\1', file) == 'Event') {
                   return(NA)
               }
               result <- gsub(pattern, '\\2', file)
               result <- strsplit(result, '\\.')[[1]][3]
               result
           },
           'channel' = {
               result <- gsub(pattern, '\\3', file)
               result <- gsub('CH', '', result)
               result
           },
           'time' = {
               result <- gsub(pattern, '\\4', file)
               milli <- gsub('(.*)_([0-9]{3})$', '\\2', result)
               result <- gsub('(.*)_([0-9]{3})$', '\\1', result)
               result <- gsub('_', '', result)
               result <- as.POSIXct(result, format='%Y%m%d%H%M%S', tz='UTC')
               result <- result + as.numeric(milli)/1e3
               result
           },
           'UTC' = {
               result <- gsub(pattern, '\\4', file)
               milli <- gsub('(.*)_([0-9]{3})$', '\\2', result)
               result <- gsub('(.*)_([0-9]{3})$', '\\1', result)
               result <- gsub('_', '', result)
               result <- as.POSIXct(result, format='%Y%m%d%H%M%S', tz='UTC')
               result <- result + as.numeric(milli)/1e3
               result
           }
    )
}
