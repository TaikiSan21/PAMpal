#' @title Plot Graphical Representations of Waveforms
#'
#' @description Fetches matching binary data from a single or multiple
#'   detections in an \linkS4class{AcousticEvent} object, then plot
#'   the resulting data
#'
#' @param x a \linkS4class{AcousticStudy} object, a list of \linkS4class{AcousticEvent}
#'   objects, or a single \linkS4class{AcousticEvent} object
#' @param UID the UID(s) of the individual detections to fetch the binary
#'   data for
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
#' @importFrom signal specgram
#' @importFrom PAMmisc wignerTransform
#' @importFrom graphics plot title
#' @export
#'
plotWaveform <- function(x, UID) {
    data <- getBinaryData(x, UID)
    for(i in seq_along(data)) {
        if(!('wave' %in% names(data[[i]]))) {
            warning('No waveform data found for UID ', names(data)[i])
            next
        }
        wav <- data[[i]]$wave
        for(c in 1:ncol(wav)) {
            plot(wav[, c], type='l')
            title(main = paste0('UID ', names(data)[i], ', Channel ', c))
        }
    }
}

#' @export
#' @rdname plotWaveform
#'
plotSpectrogram <- function(x, UID, sr=NULL, ...) {
    data <- getBinaryData(x, UID)
    for(i in seq_along(data)) {
        if(!('wave' %in% names(data[[i]]))) {
            warning('No waveform data found for UID ', names(data)[i])
            next
        }
        if(is.null(sr)) {
            if(!('sr' %in% names(data[[i]]))) {
                sr <- as.numeric(
                    readline(prompt='Sample rate not found in data, what SR should we use?')
                )
            } else {
                sr <- data[[i]]$sr
            }
        }
        wav <- data[[i]]$wave
        for(c in 1:ncol(wav)) {
            print(specgram(wav[, c], Fs = sr, ...))
            title(main = paste0('UID ', names(data)[i], ', Channel ', c))
        }
    }
}

#' @export
#' @rdname plotWaveform
#'
plotWigner <- function(x, UID, sr=NULL, ...) {
    data <- getBinaryData(x, UID)
    for(i in seq_along(data)) {
        if(!('wave' %in% names(data[[i]]))) {
            warning('No waveform data found for UID ', names(data)[i])
            next
        }
        if(is.null(sr)) {
            if(!('sr' %in% names(data[[i]]))) {
                sr <- as.numeric(
                    readline(prompt='Sample rate not found in data, what SR should we use?')
                )
            } else {
                sr <- data[[i]]$sr
            }
        }
        wav <- data[[i]]$wave
        for(c in 1:ncol(wav)) {
            wt <- wignerTransform(wav[, c], sr=sr, plot = TRUE)
            title(main = paste0('UID ', names(data)[i], ', Channel ', c))
        }
    }
}
