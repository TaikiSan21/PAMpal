#' @title Filter an AcousticStudy or AcousticEvent Object
#'
#' @description Apply dplyr-like filtering to the detecitons of an
#'   AcousticStudy or AcousticEvent object, with a special case for
#'   filtering by species for an AcousticStudy
#'
#' @param .data AcousticStudy or AcousticEvent to filter
#' @param \dots Logical expressions, syntax is identical to \link[dplyr]{filter}.
#'   There is a special case if \code{.data} is an AcousticStudy object where a
#'   logical expression using \code{species} or \code{Species} will filter by the
#'   species present in the \code{$id} of the \code{species} slot within each
#'   AcousticEvent
#' @param .preserve not used
#'
#' @return The original \code{.data} object, filtered by the given logical expressions
#'
#' @author Taiki Sakai \email{taiki.sakai@@noaa.gov}
#' @rdname filter
#' @importMethodsFrom dplyr filter
#' @importFrom rlang as_label quos parse_expr
#' @export
#'
filter.AcousticStudy <- function(.data, ..., .preserve=FALSE) {
    dotChars <- sapply(quos(...), as_label)
    isSpecies <- grepl('species|Species', dotChars)
    if(any(isSpecies)) {
        # do species filtering first
        spKeep <- rep(FALSE, length(events(.data)))
        exprText <- gsub('(species|Species)', 'species(x)$id', dotChars[isSpecies])
        for(s in seq_along(exprText)) {
            thisKeep <- sapply(events(.data), function(x) eval(parse_expr(exprText[s])))
            # browser()
            spKeep <- spKeep | thisKeep
        }
        events(.data) <- events(.data)[spKeep]
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
