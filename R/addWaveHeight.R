#' @title Add Wave Height Data to an AcousticStudy
#'
#' @description Add wave height to an AcousticStudy or AcousticEvent
#'
#' @param x an \linkS4class{AcousticStudy} to add height data to
#' @param height either a single numeric value, or a dataframe
#'   with column \code{UTC} and either column \code{waveHeight} 
#'   specifying height (m) at that time, or \code{beaufort} specifying
#'   the beaufort sea state at that time
#' @param thresh maximum time apart in seconds for matching height to
#'   data, if the closest value is more than \code{thresh} apart then the
#'   height value will be set to \code{NA}
#'
#' @details height values will be matched to the data
#'   by using data.table's rolling join with \code{roll='nearest'}. After the
#'   join is done, the time difference between the matched rows is checked
#'   and any that are greater than the set threshold are set to NA. This is
#'   done to prevent accidentally matching weird things if an incomplete set
#'   of height data is provided.
#'
#' @return the same data as \code{x}, with wave height data added. All AcousticEvents will
#'   have height data added to all detector dataframes as column \code{waveHeight}
#'
#' @author Taiki Sakai \email{taiki.sakai@@noaa.gov}
#'
#' @examples
#'
#' data(exStudy)
#' # need to update database file to local directory
#' exStudy <- addWaveHeight(exStudy, height=.5)
#' getClickData(exStudy[1])
#'
#' @importFrom data.table data.table setkeyv key setDT setDF as.data.table
#' @export
#' 
addWaveHeight <- function(x, height, thresh=3600) {
    if(is.AcousticStudy(x)) {
        for(e in seq_along(events(x))) {
            x[[e]] <- addWaveHeight(x[[e]], height)
        }
        return(x)
    }
    if(is.AcousticEvent(x)) {
        for(d in seq_along(detectors(x))) {
            thisType <- attr(x[[d]], 'calltype')
            x[[d]] <- addWaveHeight(x[[d]], height)
            attr(x[[d]], 'calltype') <- thisType
        }
        return(x)
    }
    if(is.data.frame(x)) {
        if(is.numeric(height)) {
            x$waveHeight <- height
            return(x)
        }
        if(!('UTC' %in% colnames(height))) {
            stop('Wave height data must have column "UTC"')
        }
        if(!inherits(height$UTC, 'POSIXct')) {
            height$UTC <- parseUTC(height$UTC)
        }
        if(sum(c('waveHeight', 'beaufort') %in% colnames(height)) < 1) {
            stop('Wave height data must have column "waveHeight" or "beaufort"')
        }
        setDT(height)
        if(!('waveHeight' %in% colnames(height))) {
            setkeyv(height, 'beaufort')
            height <- ppVars()$bftHeight[height, roll=-Inf]
        }
        height$heightTime <- height$UTC
        setkeyv(height, 'UTC')
        setDT(x)
        x$dataTime <- x$UTC
        setkeyv(x, 'UTC')
        x <- height[x, roll='nearest']
        x[abs(dataTime - heightTime) > thresh, c('waveHeight') := NA]
        x$UTC <- x$dataTime
        x[, c('dataTime', 'heightTime') := NULL]
        if(any(is.na(x$waveHeight))) {
            pamWarning('Some height matches exceeded time threshold, setting',
                       'value to NA.')
        }
        setDF(x)
    }
    x
}

globalVariables(c('heightTime'))
