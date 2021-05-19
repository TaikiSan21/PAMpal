#' @title Plot Graphical Representations of Waveforms
#'
#' @description Fetches matching binary data from a single or multiple
#'   detections in an \linkS4class{AcousticStudy}  object, then plot
#'   the resulting data
#'
#' @param x a \linkS4class{AcousticStudy} object, a list of \linkS4class{AcousticEvent}
#'   objects, or a single \linkS4class{AcousticEvent} object
#' @param UID the UID(s) of the individual detections to fetch the binary
#'   data for
#' @param length length of the waveform to use for plotting, in samples. The clip
#'   used will be centered around the maximum value of the waveform, if \code{length}
#'   is \code{NULL} (default), the entire waveform will be used. If \code{length} is
#'   greater than the stored clip, the waveform will be zero-padded to \code{length}
#' @param sr if \code{NULL} (default) will try to read sample rate from your
#'   data. If provided as a value will override sample rate in the data.
#' @param \dots other arguments to pass to the spectrogram or wigner functions
#'
#' @details The \code{plotSpectrogram} function uses the function
#'   \code{\link[signal]{specgram}} to plot the spectrogram, see this function
#'   for plotting options. The \code{plotWigner} function uses the function
#'   \code{\link[PAMmisc]{wignerTransform}} to plot the Wigner-Ville transform,
#'   see this function for options.
#'
#' @return Nothing, just shows plots for every channel of the waveform for
#'   each UID provided
#'
#' @author Taiki Sakai \email{taiki.sakai@@noaa.gov}
#'
#' @examples
#'
#' data(exStudy)
#' plotWaveform(exStudy, 8000003)
#' plotSpectrogram(exStudy, 8000003)
#' plotWigner(exStudy, 8000003)
#'
#' @importFrom signal specgram
#' @importFrom PAMmisc wignerTransform
#' @importFrom graphics plot title
#' @export
#'
plotWaveform <- function(x, UID, length=NULL, sr=NULL) {
    data <- getBinaryData(x,  UID, type='click')
    if(is.null(data) ||
       length(data) == 0) {
        warning('No data found for provided UID(s).')
        return(NULL)
    }
    data <- unique(data)
    for(i in seq_along(data)) {
        if(!('wave' %in% names(data[[i]]))) {
            warning('No waveform data found for UID ', names(data)[i])
            next
        }
        if(is.null(sr)) {
            if(!is.null(getSr(x))) {
                sr <- getSr(x)
            } else if(!('sr' %in% names(data[[i]]))) {
                sr <- as.numeric(
                    readline(prompt='Sample rate not found in data, what SR should we use?')
                )
            } else {
                sr <- data[[i]]$sr
            }
        }
        wav <- data[[i]]$wave
        for(c in 1:ncol(wav)) {
            if(is.null(length)) {
                plotWav <- wav[, c]
            } else {
                plotWav <- clipAroundPeak(wav[, c], length)
            }

            plot(y=plotWav, x=(1:length(plotWav))/sr*1e3, type='l', xlab = 'Time (ms)', ylab='Amplitude')
            title(main = paste0('UID ', data[[i]]$UID, ', Channel ', c))
        }
    }
}

#' @export
#' @rdname plotWaveform
#'
plotSpectrogram <- function(x, UID, length=NULL, sr=NULL, ...) {
    data <- getBinaryData(x, UID, type='click')
    if(is.null(data) ||
       length(data) == 0) {
        warning('No data found for provided UID(s).')
        return(NULL)
    }
    data <- unique(data)
    for(i in seq_along(data)) {
        if(!('wave' %in% names(data[[i]]))) {
            warning('No waveform data found for UID ', names(data)[i])
            next
        }
        if(is.null(sr)) {
            if(!is.null(getSr(x))) {
                sr <- getSr(x)
            } else if(!('sr' %in% names(data[[i]]))) {
                sr <- as.numeric(
                    readline(prompt='Sample rate not found in data, what SR should we use?')
                )
            } else {
                sr <- data[[i]]$sr
            }
        }
        wav <- data[[i]]$wave
        for(c in 1:ncol(wav)) {
            if(is.null(length)) {
                plotWav <- wav[, c]
            } else {
                plotWav <- clipAroundPeak(wav[, c], length)
            }
            print(specgram(plotWav, Fs = sr, ...))
            title(main = paste0('UID ', data[[i]]$UID, ', Channel ', c))
        }
    }
}

#' @export
#' @rdname plotWaveform
#'
plotWigner <- function(x, UID, length=NULL, sr=NULL, ...) {
    data <- getBinaryData(x, UID, type='click')
    if(is.null(data) ||
       length(data) == 0) {
        warning('No data found for provided UID(s).')
        return(NULL)
    }
    data <- unique(data)
    for(i in seq_along(data)) {
        if(!('wave' %in% names(data[[i]]))) {
            warning('No waveform data found for UID ', names(data)[i])
            next
        }
        if(is.null(sr)) {
            if(!is.null(getSr(x))) {
                sr <- getSr(x)
            } else if(!('sr' %in% names(data[[i]]))) {
                sr <- as.numeric(
                    readline(prompt='Sample rate not found in data, what SR should we use?')
                )
            } else {
                sr <- data[[i]]$sr
            }
        }
        wav <- data[[i]]$wave
        for(c in 1:ncol(wav)) {
            if(is.null(length)) {
                plotWav <- wav[, c]
            } else {
                plotWav <- clipAroundPeak(wav[, c], length)
            }
            wt <- wignerTransform(plotWav, sr=sr, plot = TRUE)
            title(main = paste0('UID ', data[[i]]$UID, ', Channel ', c))
        }
    }
}
