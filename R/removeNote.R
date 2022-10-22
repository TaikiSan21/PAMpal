#' @title removeNote
#'
#' @description Remove a note added with \link{addNote}
#'
#' @param x An \linkS4class{AcousticStudy} or \linkS4class{AcousticEvent} object
#' @param index The index of the note to remove, order matches the output of \link{getNotes}
#'
#' @return For \code{addNote}, the same data as \code{x}, with notes added.
#'   For \code{getNotes}, a list of all notes present in \code{x}
#'
#' @author Taiki Sakai \email{taiki.sakai@@noaa.gov}
#'
#' @examples
#'
#' data(exStudy)
#' exStudy <- addNote(exStudy, to='study', label='Note1',
#'                    note='My first note for this study')
#' exStudy <- addNote(exStudy, to='event', evNum=1:2, label='Note2',
#'                    note='A note for the first two events')
#' exStudy
#' removeNote(exStudy, 1)
#' removeNote(exStudy, 2)
#' removeNote(exStudy, 3)
#'
#' @export
#'
removeNote <- function(x, index) {
    if(index > length(unlist(getNotes(x)))) {
        return(x)
    }
    nStudy <- length(ancillary(x)$notes)
    if(is.AcousticEvent(x) ||
       index <= nStudy) {
        ancillary(x)$notes[index] <- NULL
        if(length(ancillary(x)$notes) == 0) {
            ancillary(x)$notes <- NULL
        }
        return(x)
    }
    index <- index - nStudy
    nEvent <- sapply(events(x), function(e) length(ancillary(e)$notes))
    cumEv <- cumsum(nEvent)
    whichEv <- min(which(cumEv >= index))
    noteIx <- ifelse(whichEv == 1, index, index - cumEv[whichEv-1])
    ancillary(x[[whichEv]])$notes[noteIx] <- NULL
    if(length(ancillary(x[[whichEv]])$notes) == 0) {
        ancillary(x[[whichEv]])$notes <- NULL
    }
    x
}

