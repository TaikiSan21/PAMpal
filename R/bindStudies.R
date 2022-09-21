#' @title Combine AcousticStudy Objects
#'
#' @details All events will be combined into one large list of events. Files,
#'   settings, effort, models, GPS, and ancillary fields will be combined
#'   using the \link[PAMmisc]{squishList} function from the PAMmisc package
#'   (dataframes are combined, vectors are appended). The id is changed by
#'   pasting all IDs together along with a note that they have been combined.
#'   Note that the \linkS4class{PAMpalSettings} object in the pps slot is just
#'   left as the \code{pps} in the first AcousticStudy to be combined, and thus
#'   is not representative of the new combined AcousticStudy
#'
#' @description Combines multiple AcousticStudy objects (or lists of these)
#'   into a single object
#'
#' @param \dots AcousticStudy objects, or a list of AcousticStudy objects
#'
#' @return A single AcousticStudy object
#'
#' @author Taiki Sakai \email{taiki.sakai@@noaa.gov}
#'
#' @importFrom rlang list2
#' @importFrom purrr discard
#'
#' @export
#'
bindStudies <- function(...) {
    dots <- rlang::list2(...)
    if(is.list(dots)) {
        dots <- unlist(dots, recursive = FALSE)
    }
    dots <- lapply(dots, function(x) {
        if(is.AcousticStudy(x)) {
            return(x)
        }
        if(is.list(x)) {
            return(unlist(x, recursive = FALSE))
        }
        NULL
    })
    dots <- purrr::discard(dots, is.null)
    if(!all(sapply(dots, is.AcousticStudy))) {
        stop('Inputs must all be AcousticStudy objects or lists of AcousticStudy objects.')
    }
    newAc <- dots[[1]]
    newEv <- unlist(lapply(dots, events), recursive = FALSE)
    evNames <- unique(names(newEv))
    if(length(evNames) < length(newEv)) {
        newNames <- c()
        for(e in seq_along(evNames)) {
            whichThis <- which(names(newEv) == evNames[e])
            if(length(whichThis) == 1) {
                next
            }
            for(w in seq_along(whichThis)) {
                id(newEv[[whichThis[w]]]) <- paste0(names(newEv)[whichThis[w]], '_', w)
                names(newEv)[whichThis[w]] <- id(newEv[[whichThis[w]]])
                newNames <- c(newNames, id(newEv[[whichThis[w]]]))
            }
        }
        pamWarning('Duplicate names found in combined study, numbers have been added to the end ',
                   'to create new names. New names are ', newNames)
    }

    events(newAc) <- newEv
    newGps <- doAcCombine(dots, gps, unique=TRUE, df=TRUE)
    if(is.list(newGps) && length(newGps) == 0) {
        newGps <- data.frame()
    }
    gps(newAc) <- newGps
    files(newAc) <- doAcCombine(dots, files, unique=TRUE)
    settings(newAc) <- doAcCombine(dots, settings, unique=TRUE)
    newEffort <- doAcCombine(dots, effort, df=TRUE)
    if(is.list(newEffort) && length(newEffort) == 0) {
        newEffort <- data.frame()
    }
    effort(newAc) <- newEffort
    models(newAc) <- doAcCombine(dots, models)
    ancillary(newAc) <- doAcCombine(dots, ancillary)
    ancillary(newAc)$processDate <- as.POSIXct(ancillary(newAc)$processDate, origin='1970-01-01 00:00:00', tz='UTC')
    id(newAc) <- paste0('COMBINED: ', paste0(sapply(dots, id), collapse=', '))
    newAc <- .addPamWarning(newAc)
    newAc
}

#' @importFrom PAMmisc squishList
#'
doAcCombine <- function(dots, FUN, unique=FALSE, df=TRUE) {
    data <- lapply(dots, FUN)
    if(df) {
        return(bind_rows(data))
    }
    squishList(data, unique=unique)
}
