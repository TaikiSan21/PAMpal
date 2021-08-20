#' @title Filter an AcousticStudy or AcousticEvent Object
#'
#' @description Apply dplyr-like filtering to the detecitons of an
#'   AcousticStudy or AcousticEvent object, with a special case for
#'   filtering by species for an AcousticStudy
#'
#' @param .data \linkS4class{AcousticStudy} or \linkS4class{AcousticEvent} to filter
#' @param \dots Logical expressions, syntax is identical to \link[dplyr]{filter}.
#'   There are special cases to filter by environmental variables, species ID, 
#'   database, or detector name. See details.
#' @param .preserve not used
#'
#' @details Most expression provided will be used to filter out detections based on
#'   calculated parameters. 
#'   
#'   If the name of an environmental variable added using
#'   \link{matchEnvData} is provided, will filter to only events with environmental
#'   variables matching those conditions. 
#'   
#'   If a provided logical expression uses
#'   \code{"species"} or \code{"Species"}, then events will be filtered using the
#'   species present in the \code{$id} of the \code{species} slot of each event.
#'   
#'   If a provided logical expression uses \code{"database"} or \code{"Database"}, 
#'   then only events with databases matching the expression in \code{files(.data)$db}
#'   will remain
#'   
#'   If a provided logical expression uses \code{"detector"} or \code{"Detector"}, then
#'   only detections from detectors with names matching the expression will remain in
#'   events. Any events left with no detections will be removed.
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
        pamWarning('Did you put "=" when you meant "=="? This filter will not be applied.')
    }
    # do event level filters first
    # browser()
    # isSpecies <- grepl('^species|^Species', dotChars)
    isSpecies <- grepl('species|Species', dotChars)
    if(any(isSpecies)) {
        # do species filtering first
        naSp <- sapply(events(.data), function(x) is.na(species(x)$id))
        if(any(naSp)) {
            pamWarning('Attempting to filter by species, but ', sum(naSp),
                    ' species have not been set. These will be removed from',
                    ' the filtered results.')
        }
        spKeep <- rep(TRUE, length(events(.data)))
        # exprText <- gsub('(^species|^Species)', 'species(x)$id', dotChars[isSpecies])
        exprText <- gsub('^(.*?)species(.*)', '\\1species(x)$id\\2', dotChars[isSpecies], ignore.case=TRUE)
        for(s in seq_along(exprText)) {
            thisKeep <- sapply(events(.data), function(x) eval(parse_expr(exprText[s])))
            thisKeep[is.na(thisKeep)] <- FALSE
            # browser()
            spKeep <- spKeep & thisKeep
        }
        events(.data) <- events(.data)[spKeep]
        if(length(events(.data)) == 0) {
            .data <- .addPamWarning(.data)
            return(.data)
        }
    }

    isDb <- grepl('database|Database', dotChars)
    if(any(isDb)) {
        dbKeep <- rep(TRUE, length(events(.data)))
        exprText <- gsub('^(.*?)database(.*)', '\\1files(x)$db\\2', dotChars[isDb], ignore.case=TRUE)
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
            .data <- .addPamWarning(.data)
            return(.data)
        }
    }
    # do enviro?
    if(!is.null(ancillary(.data[[1]])$environmental)) {
        envNames <- names(ancillary(.data[[1]])$environmental)
        envNames <- envNames[!(envNames %in% c('UTC', 'Longitude', 'Latitude'))]
        hasEnv <- sapply(dotChars, function(d) any(sapply(envNames, function(nm) grepl(nm, d))))
        if(any(hasEnv)) {
            evDf <- bind_rows(lapply(events(.data), function(x) {
                ancillary(x)$environmental
            }))

            filteredEv <- doFilter(evDf[, !(names(evDf) %in% c('UTC', 'Longitude', 'Latitude')), drop=FALSE], ...)
            events(.data) <- events(.data)[names(events(.data)) %in% unique(filteredEv$event)]
            if(length(events(.data)) == 0) {
                .data <- .addPamWarning(.data)
                return(.data)
            }
        }
    }

    events(.data) <- lapply(events(.data), function(x) {
        filter(x, ...)
    })
    isNull <- sapply(events(.data), is.null)
    events(.data) <- events(.data)[!isNull]
    .data <- .addPamWarning(.data)
    .data
}

#' @export
#'
filter.AcousticEvent <- function(.data, ..., .preserve=FALSE) {
    # browser()
    dotChars <- sapply(quos(...), as_label)
    # isDetector <- grepl('^.{0,3}detector|^.{0,3}Detector', dotChars)
    isDetector <- grepl('\\b[Dd]etector\\b', dotChars, ignore.case=TRUE)
    detKeep <- rep(TRUE, length(detectors(.data)))
    if(any(isDetector)) {
        exprText <- gsub('^(.*?)\\bdetector\\b(.*)', '\\1names(detectors(.data))\\2', dotChars[isDetector], ignore.case=TRUE)
        for(s in seq_along(exprText)) {
            thisKeep <- eval(parse_expr(exprText[s]))
            thisKeep[is.na(thisKeep)] <- FALSE
            detKeep <- detKeep & thisKeep
        }
        if(!any(detKeep)) {
            return(NULL)
        }
        detectors(.data) <- detectors(.data)[detKeep]
    }
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
