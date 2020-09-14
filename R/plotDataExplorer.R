#' @title Explore Data in an Interactive Plot
#'
#' @description Creates an interactive plot of detector data. Allows user to
#'   choose which numeric data to plot, and will allow user to both color and
#'   facet the plot by any columns that are characters or factors
#'
#' @param x data to plot, can be an \code{AcousticStudy}, \code{AcousticEvent},
#'   data.frame or a list of \code{AcousticEvent} objects
#' @param maxCategories maximum number of categories to color and facet by. Only
#'   character and factor data with a number of unique values less than or equal
#'   to this number will be shown as options for selecting colors and facets. Not
#'   recommended to increase this value much beyond 20, trying to plot a large number
#'   of colors will cause R to be sad.
#' @param callType the specific type of call to plot. If \code{NULL} (default),
#'   will prompt user to choose which type if more than one is present.
#'
#' @return nothing, just plots
#'
#' @author Taiki Sakai \email{taiki.sakai@@noaa.gov}
#'
#' @importFrom manipulate manipulate
#' @importFrom ggplot2 ggplot geom_density aes_string facet_wrap
#' @export
#'
plotDataExplorer <- function(x, maxCategories=15, callType=NULL) {
    globalVariables(c('dataPicker', 'colPicker', 'fCheck'))
    x <- getDetectorData(x)
    if(is.list(x)) {
        if(length(x) == 1) {
            x <- x[[1]]
        } else if(length(x) > 1) {
            if(!is.null(callType) &&
               callType %in% names(x)) {
                x <- x[[callType]]
            } else {
                choice <- menu(choices = names(x), title = 'What call type would you like to explore?')
                x <- x[[choice]]
            }
        } else { # case length 0 data, shouldn't happen
            stop('No data found')
        }
    }
    if(!is.data.frame(x)) { # this shouldnt be possible
        stop('Could not convert "x" to a data frame.')
    }
    charCols <- colnames(x)[sapply(x, function(c) is.character(c) || is.factor(c))]
    nColors <- sapply(charCols, function(c) {
        length(unique(x[[c]]))
    })
    charCols <- charCols[nColors <= maxCategories]
    numCols <- colnames(x)[sapply(x, function(c) is.numeric(c) || inherits(c, 'POSIXct'))]

    dPick <- myPicker(numCols, label = 'Data')
    cPick <- myPicker(charCols, none=TRUE, label = 'Color')
    facetPick <- myPicker(charCols, none=TRUE, label = 'Facet By')

    manipulate({
        ggplot(x) +
            geom_density(aes_string(x=dataPicker, col=colPicker)) +
            facet_wrap(facets = fCheck)

    },
    dataPicker =dPick,
    colPicker = cPick,
    fCheck = facetPick)
}

globalVariables(c('dataPicker', 'colPicker', 'fCheck'))

myPicker <- function(x, none=FALSE, label) {
    vals <- as.list(x)
    if(none) {
        vals <- c('None' = list(NULL), vals)
        x <- c('None', x)
    }
    names(vals) <- x

    myPick <- list(type = 1,
                   choices = x,
                   values = vals,
                   initialValue = x[1],
                   label = label)
    class(myPick) <- 'manipulator.picker'
    myPick
}
