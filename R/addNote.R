#' @title addNote
#'
#' @description Adds a note to an AcousticEvent or AcousticStudy. Notes can either
#'   be accessed with the "getNotes" function, or up to 6 notes will be printed
#'   when the object is printed
#'
#' @param x An \linkS4class{AcousticStudy} or \linkS4class{AcousticEvent} object
#' @param to One of "study" or "event", which object to add the note to
#' @param evNum If \code{x} is an AcousticStudy and \code{to} is "event", the
#'   number or name of the event(s) to add notes to (can be a vector of numbers
#'   or names to add the same note to multiple events)
#' @param label (optional) a short header or label for the note. Recommended to
#'   set this as a sumamry of the more detailed note
#' @param note the full note message
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
#' exStudy <- addNote(exStudy[[1]], to='event', label='Note3',
#'                    note='A second note for the first event')
#' exStudy
#'
#' @export
#'
addNote <- function(x, to=c('study', 'event'), evNum=1, label=NULL, note) {
    to <- match.arg(to)
    isStudy <- is.AcousticStudy(x)
    isEvent <- is.AcousticEvent(x)
    if(!any(c(isStudy, isEvent))) {
        stop('x must be an AcousticEvent or AcousticStudy')
    }
    if(to == 'study' &&
       isEvent) {
        stop('"x" is an AcousticEvent, but "to" is set to "study"')
    }
    if(to == 'event' &&
       isStudy) {
        for(e in evNum) {
            x[[e]] <- doNoteAdd(x[[e]], label, note)
        }
        return(x)
    }
    x <- doNoteAdd(x, label, note)
    x
}

doNoteAdd <- function(x, label, note) {
    if(is.null(label)) {
        ancillary(x)$notes <- c(ancillary(x)$notes, note)
    } else {
        ancillary(x)$notes[[label]] <- note
    }
    x
}

formatNotes <- function(x, nchar=70, n=6, nSpace=0) {
    if(is.AcousticEvent(x) ||
       is.AcousticStudy(x)) {
        x <- getNotes(x)
    }
    if(is.character(x)) {
        if(nchar(x) > nchar) {
            x <- paste0(substr(x, 1, nchar), '...')
        }
        return(paste0(paste0(rep(' ', nSpace), collapse=''), x))
    }
    nx <- length(unlist(x))
    end <- ''
    if(nx > n) {
        x <- trimToN(x, n)
        end <- paste0('\n(', nx-n, ' more note(s) not shown)')
    }
    notes <- vector('list', length=length(x))
    for(i in seq_along(x)) {
        noteName <- names(x[i])
        note <- formatNotes(x[[i]], nchar, n, nSpace=nSpace+2)

        if(noteName != '') {
            noteName <- paste0(paste0(rep(' ', nSpace), collapse=''), noteName, ':\n')
        } else {
            note <- substr(note, 3, nchar+nSpace+2+3)
        }
        note[1] <- paste0(noteName, note[1])
        # note[length(note)] <- paste0(note[length(note)], end)
        notes[[i]] <- note
    }
    # printN(notes, collapse='\n')
    notes <- unlist(notes)
    notes[length(notes)] <- paste0(notes[length(notes)], end)
    notes
}
#' @rdname addNote
#' @export
#'
getNotes <- function(x) {
    if(is.AcousticEvent(x)) {
        out <- ancillary(x)$notes
        if(length(out) == 0) {
            return(NULL)
        }
        return(out)
    }
    if(is.AcousticStudy(x)) {
        out <- list(studyNotes = ancillary(x)$notes)
        evNotes <- lapply(events(x), getNotes)
        evNotes <- evNotes[sapply(evNotes, function(e) !is.null(e))]
        if(length(evNotes) > 0) {
            out$eventNotes <- evNotes
        }
        if(length(out$studyNotes) == 0) {
            out$studyNotes <- NULL
        }
        if(is.null(out$studyNotes) &&
           length(evNotes) == 0) {
            return(NULL)
        }
    }
    out
}
# auto add note to study if getWarnings > 0
### PROBLEM IS THAT THERE IS NO POINT WHERE WE CHECK IF WERE IN
# IN A VECTOR AND JUST CUT IT OFF TO N AT THAT POINT STILL FKIN WEIRD
trimToN <- function(x, n=6) {
    if(n < 1) {
        return(NULL)
    }
    if(!is.list(x[[1]])) {
        x <- x[1:n]
        return(x)
    }
    if(length(unlist(x)) <= n) {
        return(x)
    }
    if(length(unlist(x[[1]])) > n) {
        x[[1]] <- trimToN(x[[1]], n)
        return(x[1])
    }
    x[2:length(x)] <- trimToN(x[2:length(x)], n - length(unlist(x[[1]])))
    x
}
