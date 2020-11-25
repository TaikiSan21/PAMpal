#' @title Filter an AcousticStudy or AcousticEvent Object
#'
#' @description Apply dplyr-like filtering to the detecitons of an
#'   AcousticStudy or AcousticEvent object, with a special case for
#'   filtering by species for an AcousticStudy
#'
#' @param .data \linkS4class{AcousticStudy} or \linkS4class{AcousticEvent} to filter
#' @param \dots Logical expressions, syntax is identical to \link[dplyr]{filter}.
#'   There is a special case if \code{.data} is an AcousticStudy object where a
#'   logical expression using \code{species} or \code{Species} will filter by the
#'   species present in the \code{$id} of the \code{species} slot within each
#'   AcousticEvent
#' @param .preserve not used
#'
#' @return The original \code{.data} object, filtered by the given logical expressions
#'
#' @examples
#'
#' # create example data
#' data(exStudy)
#' exStudy <- setSpecies(exStudy, method='manual', value=letters[1:2])
#' filterData <- filter(exStudy, peak < 20)
#' getDetectorData(filterData)$click
#'
#' filterData <- filter(exStudy, species == 'a')
#' species(filterData[[1]])
#'
#' @author Taiki Sakai \email{taiki.sakai@@noaa.gov}
#' @rdname filter
#' @importFrom dplyr filter
#' @importFrom rlang as_label quos parse_expr
#' @export
#'
filter.AcousticStudy <- function(.data, ..., .preserve=FALSE) {
    dotChars <- sapply(quos(...), as_label)
    notFilt <- names(dotChars) != ''
    if(any(notFilt)) {
        warning('Did you put "=" when you meant "=="? This filter will not be applied.')
    }
    # do event level filters first
    isSpecies <- grepl('^species|^Species', dotChars)
    if(any(isSpecies)) {
        # do species filtering first
        naSp <- sapply(events(.data), function(x) is.na(species(x)$id))
        if(any(naSp)) {
            warning('Attempting to filter by species, but ', sum(naSp),
                    ' species have not been set. These will be removed from',
                    ' the filtered results.')
        }
        spKeep <- rep(TRUE, length(events(.data)))
        exprText <- gsub('(^species|^Species)', 'species(x)$id', dotChars[isSpecies])
        for(s in seq_along(exprText)) {
            thisKeep <- sapply(events(.data), function(x) eval(parse_expr(exprText[s])))
            thisKeep[is.na(thisKeep)] <- FALSE
            # browser()
            spKeep <- spKeep & thisKeep
        }
        events(.data) <- events(.data)[spKeep]
        if(length(events(.data)) == 0) {
            return(.data)
        }
    }

    isDb <- grepl('^database|^Database', dotChars)
    if(any(isDb)) {
        dbKeep <- rep(TRUE, length(events(.data)))
        exprText <- gsub('(^database|^Database)', 'files(x)$db', dotChars[isDb])
        studyExpr <- gsub('\\(x\\)', '\\(\\.data\\)', exprText)
        for(d in seq_along(exprText)) {
            thisKeep <- sapply(events(.data), function(x) eval(parse_expr(exprText[d])))
            dbKeep <- dbKeep & thisKeep
        }
        studyKeep <- rep(TRUE, length(files(.data)$db))
        for(s in studyExpr) {
            thisKeep <- eval(parse_expr(s))
            if(all(is.logical(thisKeep))) {
                studyKeep <- studyKeep & thisKeep
            }
        }
        events(.data) <- events(.data)[dbKeep]
        files(.data)$db <- files(.data)$db[studyKeep]
        if(length(events(.data)) == 0) {
            return(.data)
        }
    }
    # do enviro?
    if(!is.null(ancillary(.data[[1]])$environmental)) {
        evDf <- bind_rows(lapply(events(.data), function(x) {
            ancillary(x)$environmental
        }))
        filteredEv <- doFilter(evDf[, !(names(evDf) %in% c('UTC', 'Longitude', 'Latitude')), drop=FALSE], ...)
        events(.data) <- events(.data)[unique(filteredEv$event)]
        if(length(events(.data)) == 0) {
            return(.data)
        }
    }

    events(.data) <- lapply(events(.data), function(x) {
        filter(x, ...)
    })
    isNull <- sapply(events(.data), is.null)
    events(.data) <- events(.data)[!isNull]
    .data
}

#' @export
#'
filter.AcousticEvent <- function(.data, ..., .preserve=FALSE) {
    detectors(.data) <- lapply(detectors(.data), function(x) {
        doFilter(x, ...)
    })
    detNums <- sapply(detectors(.data), nrow)
    if(all(detNums == 0)) {
        return(NULL)
    }
    detectors(.data) <- detectors(.data)[detNums > 0]
    .data
}

doFilter <- function(.x, ...) {
    dotChars <- sapply(quos(...), as_label)
    hasCol <- sapply(dotChars, function(d) any(sapply(colnames(.x), function(c) grepl(c, d))))
    if(!any(hasCol)) {
        return(.x)
    }
    filter(.x, !!!quos(...)[hasCol])
}
