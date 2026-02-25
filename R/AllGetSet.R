# All get/set
# Get/Set for AcousticEvent class -----------------------------------------

#' @title \code{AcousticEvent} and \code{AcousticStudy} accessors
#'
#' @description Accessors for slots in \linkS4class{AcousticEvent}
#'   and \linkS4class{AcousticStudy} objects
#'
#' @param x a \linkS4class{AcousticEvent} or \linkS4class{AcousticStudy} object
#' @param value value to assign with accessor
#' @param i index of the object to access
#' @param j not used
#' @param name name of the object to access
#' @param \dots other arguments to pass to methods
#'
#' @return
#' \describe{
#'   \item{id}{a unique id or name for this object}
#'   \item{settings}{a named list of settings for each detector and localization or recorder}
#'   \item{detectors}{a list of detector data frames}
#'   \item{localizations}{list of localizations}
#'   \item{species}{list of species classifications}
#'   \item{files}{list of files used to create this object}
#'   \item{events}{a list of \linkS4class{AcousticEvent} objects}
#'   \item{gps}{a dataframe containing gps data}
#'   \item{pps}{the \linkS4class{PAMpalSettings} object used to create this}
#'   \item{effort}{something about effort?}
#'   \item{ancillary}{miscellaneous extra data}
#'   }
#'
#' @author Taiki Sakai \email{taiki.sakai@@noaa.gov}
#'
#' @name PAMpal.accessors
#'
#' @importFrom methods setGeneric setMethod validObject
#'
NULL

#' @export
#' @rdname PAMpal.accessors
#'
setGeneric('settings', function(x, ...) standardGeneric('settings'))

#' @export
#' @rdname PAMpal.accessors
#' @aliases settings
#'
setMethod('settings', 'AcousticEvent', function(x, ...) x@settings)

#' @export
#' @rdname PAMpal.accessors
#'
setGeneric('settings<-', function(x, value) standardGeneric('settings<-'))

#' @export
#' @rdname PAMpal.accessors
#' @aliases settings
#'
setMethod('settings<-', 'AcousticEvent', function(x, value) {
    x@settings <- value
    validObject(x)
    x
})

#' @export
#' @rdname PAMpal.accessors
#'
setGeneric('localizations', function(x, ...) standardGeneric('localizations'))

#' @export
#' @rdname PAMpal.accessors
#' @aliases localizations
#'
setMethod('localizations', 'AcousticEvent', function(x, ...) x@localizations)

#' @export
#' @rdname PAMpal.accessors
#'
setGeneric('localizations<-', function(x, value) standardGeneric('localizations<-'))

#' @export
#' @rdname PAMpal.accessors
#' @aliases localizations
#'
setMethod('localizations<-', 'AcousticEvent', function(x, value) {
    x@localizations <- value
    validObject(x)
    x
})

#' @export
#' @rdname PAMpal.accessors
#'
setGeneric('id', function(x, ...) standardGeneric('id'))

#' @export
#' @rdname PAMpal.accessors
#' @aliases id
#'
setMethod('id', 'AcousticEvent', function(x, ...) x@id)

#' @export
#' @rdname PAMpal.accessors
#'
setGeneric('id<-', function(x, value) standardGeneric('id<-'))

#' @export
#' @rdname PAMpal.accessors
#' @aliases id
#'
setMethod('id<-', 'AcousticEvent', function(x, value) {
    x@id <- value
    validObject(x)
    x
})

#' @export
#' @rdname PAMpal.accessors
#'
setGeneric('detectors', function(x, ...) standardGeneric('detectors'))

#' @export
#' @rdname PAMpal.accessors
#' @aliases detectors
#'
setMethod('detectors', 'AcousticEvent', function(x, ...) x@detectors)

#' @export
#' @rdname PAMpal.accessors
#'
setGeneric('detectors<-', function(x, value) standardGeneric('detectors<-'))

#' @export
#' @rdname PAMpal.accessors
#' @aliases detectors
#'
setMethod('detectors<-', 'AcousticEvent', function(x, value) {
    x@detectors <- value
    validObject(x)
    x
})

#' @export
#' @rdname PAMpal.accessors
#'
setGeneric('species', function(x, ...) standardGeneric('species'))

#' @export
#' @rdname PAMpal.accessors
#' @aliases species
#'
setMethod('species', 'AcousticEvent', function(x, ...) x@species)

#' @export
#' @param type species type to select
#' @rdname PAMpal.accessors
#' @aliases species
#'
setMethod('species', 'AcousticStudy', function(x, type='id',...) {
    sapply(events(x), function(e) species(e)[[type]])
})

#' @export
#' @rdname PAMpal.accessors
#'
setGeneric('species<-', function(x, value) standardGeneric('species<-'))

#' @export
#' @rdname PAMpal.accessors
#' @aliases species
#'
setMethod('species<-', 'AcousticEvent', function(x, value) {
    x@species <- value
    validObject(x)
    x
})

#' @export
#' @rdname PAMpal.accessors
#'
setGeneric('files', function(x, ...) standardGeneric('files'))

#' @export
#' @rdname PAMpal.accessors
#' @aliases files
#'
setMethod('files', 'AcousticEvent', function(x, ...) x@files)

#' @export
#' @rdname PAMpal.accessors
#'
setGeneric('files<-', function(x, value) standardGeneric('files<-'))

#' @export
#' @rdname PAMpal.accessors
#' @aliases files
#'
setMethod('files<-', 'AcousticEvent', function(x, value) {
    x@files <- value
    validObject(x)
    x
})

#' @export
#' @rdname PAMpal.accessors
#'
setGeneric('ancillary', function(x, ...) standardGeneric('ancillary'))

#' @export
#' @rdname PAMpal.accessors
#' @aliases ancillary
#'
setMethod('ancillary', 'AcousticEvent', function(x, ...) x@ancillary)

#' @export
#' @rdname PAMpal.accessors
#'
setGeneric('ancillary<-', function(x, value) standardGeneric('ancillary<-'))

#' @export
#' @rdname PAMpal.accessors
#' @aliases ancillary
#'
setMethod('ancillary<-', 'AcousticEvent', function(x, value) {
    x@ancillary <- safeListAdd(x@ancillary, value)
    validObject(x)
    x
})

#' @export
#' @rdname PAMpal.accessors
#'
setMethod('[', signature(x='AcousticEvent', i='ANY', j='ANY'), function(x, i) {
    x@detectors[i]
})

#' @export
#' @rdname PAMpal.accessors
#'
setMethod('[<-', 'AcousticEvent', function(x, i, value) {
    x@detectors[i] <- value
    validObject(x)
    x
})

#' @export
#' @rdname PAMpal.accessors
#'
setMethod('$', 'AcousticEvent', function(x, name) {
    '[['(x@detectors, name)
})

#' @export
#' @rdname PAMpal.accessors
#'
setMethod('$<-', 'AcousticEvent', function(x, name, value) {
    x@detectors[[name]] <- value
    validObject(x)
    x
})

#' @export
#' @rdname PAMpal.accessors
#'
setMethod('[[', signature(x='AcousticEvent', i='ANY', j='ANY'), function(x, i, ...) {
    '[['(x@detectors, i)
})

#' @export
#' @rdname PAMpal.accessors
#'
setMethod('[[<-', 'AcousticEvent', function(x, i, value) {
    x@detectors[[i]] <- value
    validObject(x)
    x
})

#' @importFrom utils .DollarNames
#' @export
#'
.DollarNames.AcousticEvent <- function(x, pattern='') {
    grep(pattern, names(detectors(x)), value=TRUE)
}

#  Get/Set for AcousticStudy class -----------------------------------------------
# #' @title \code{AcousticStudy} accessors
# #'
# #' @description Accessors for slots in \linkS4class{AcousticStudy} objects
#'
# #' @param x a \linkS4class{AcousticStudy} object
# #' @param value value to assign with accessor
# #' @param \dots other arguments to pass to methods
#'
# #' @return
# #' \describe{
# #'   \item{id}{a name or id for this study}
# #'   \item{events}{a list of \linkS4class{AcousticEvent} objects}
# #'   \item{files}{a list of files}
# #'   \item{gps}{a dataframe containing gps data}
# #'   \item{pps}{the \linkS4class{PAMpalSettings} object used to create this}
# #'   \item{settings}{a named list of settings for each detector and localization}
# #'   \item{effort}{something about effort?}
# #'   \item{ancillary}{miscellaneous extra data}
# #'   }
#'
# #' @author Taiki Sakai \email{taiki.sakai@@noaa.gov}
#'
# #' @name PAMpal.accessors
#'
# #' @importFrom methods setGeneric setMethod validObject
#'
# NULL

#' @export
#' @rdname PAMpal.accessors
#' @aliases id
#'
setMethod('id', 'AcousticStudy', function(x, ...) x@id)

#' @export
#' @rdname PAMpal.accessors
#' @aliases id
#'
setMethod('id<-', 'AcousticStudy', function(x, value) {
    x@id <- value
    validObject(x)
    x
})

#' @export
#' @rdname PAMpal.accessors
#' @aliases files
#'
setMethod('files', 'AcousticStudy', function(x, ...) x@files)

#' @export
#' @rdname PAMpal.accessors
#' @aliases files
#'
setMethod('files<-', 'AcousticStudy', function(x, value) {
    x@files <- value
    validObject(x)
    x
})

#' @export
#' @rdname PAMpal.accessors
#'
setGeneric('gps', function(x, ...) standardGeneric('gps'))

#' @export
#' @rdname PAMpal.accessors
#' @aliases gps
#'
setMethod('gps', 'AcousticStudy', function(x, ...) x@gps)

#' @export
#' @rdname PAMpal.accessors
#'
setGeneric('gps<-', function(x, value) standardGeneric('gps<-'))

#' @export
#' @rdname PAMpal.accessors
#' @aliases gps
#'
setMethod('gps<-', 'AcousticStudy', function(x, value) {
    x@gps <- value
    validObject(x)
    x
})

#' @export
#' @rdname PAMpal.accessors
#' @aliases detectors
#'
setMethod('detectors', 'AcousticStudy', function(x, ...) {
    getDetectorData(x)
})

#' @export
#' @rdname PAMpal.accessors
#'
setGeneric('events', function(x, ...) standardGeneric('events'))

#' @export
#' @rdname PAMpal.accessors
#' @aliases events
#'
setMethod('events', 'AcousticStudy', function(x, ...) x@events)

#' @export
#' @rdname PAMpal.accessors
#'
setGeneric('events<-', function(x, value) standardGeneric('events<-'))

#' @export
#' @rdname PAMpal.accessors
#' @aliases events
#'
setMethod('events<-', 'AcousticStudy', function(x, value) {
    x@events <- value
    validObject(x)
    x
})

#' @export
#' @rdname PAMpal.accessors
#' @aliases settings
#'
setMethod('settings', 'AcousticStudy', function(x, ...) x@settings)

#' @export
#' @rdname PAMpal.accessors
#' @aliases settings
#'
setMethod('settings<-', 'AcousticStudy', function(x, value) {
    x@settings <- value
    validObject(x)
    x
})

#' @export
#' @rdname PAMpal.accessors
#'
setGeneric('effort', function(x, ...) standardGeneric('effort'))

#' @export
#' @rdname PAMpal.accessors
#' @aliases effort
#'
setMethod('effort', 'AcousticStudy', function(x, ...) x@effort)

#' @export
#' @rdname PAMpal.accessors
#'
setGeneric('effort<-', function(x, value) standardGeneric('effort<-'))

#' @export
#' @rdname PAMpal.accessors
#' @aliases effort
#'
setMethod('effort<-', 'AcousticStudy', function(x, value) {
    x@effort <- value
    validObject(x)
    x
})

#' @export
#' @rdname PAMpal.accessors
#'
setGeneric('pps', function(x, ...) standardGeneric('pps'))

#' @export
#' @rdname PAMpal.accessors
#' @aliases pps
#'
setMethod('pps', 'AcousticStudy', function(x, ...) x@pps)

#' @export
#' @rdname PAMpal.accessors
#'
setGeneric('pps<-', function(x, value) standardGeneric('pps<-'))

#' @export
#' @rdname PAMpal.accessors
#' @aliases pps
#'
setMethod('pps<-', 'AcousticStudy', function(x, value) {
    x@pps <- value
    validObject(x)
    x
})

#' @export
#' @rdname PAMpal.accessors
#' @aliases ancillary
#'
setMethod('ancillary', 'AcousticStudy', function(x, ...) x@ancillary)

#' @export
#' @rdname PAMpal.accessors
#' @aliases ancillary
#'
setMethod('ancillary<-', 'AcousticStudy', function(x, value) {
    x@ancillary <- safeListAdd(x@ancillary, value)
    validObject(x)
    x
})

#' @export
#' @rdname PAMpal.accessors
#'
setGeneric('models', function(x, ...) standardGeneric('models'))

#' @export
#' @rdname PAMpal.accessors
#' @aliases models
#'
setMethod('models', 'AcousticStudy', function(x, ...) x@models)

#' @export
#' @rdname PAMpal.accessors
#'
setGeneric('models<-', function(x, value) standardGeneric('models<-'))

#' @export
#' @rdname PAMpal.accessors
#' @aliases models
#'
setMethod('models<-', 'AcousticStudy', function(x, value) {
    x@models <- value
    validObject(x)
    x
})

#' @export
#' @rdname PAMpal.accessors
#'
setMethod('[', signature(x='AcousticStudy', i='ANY', j='ANY'), function(x, i) {
    x@events <- x@events[i]
    x@events <- x@events[sapply(x@events, function(e) {
        !is.null(e)
    })]
    x
})

#' @export
#' @rdname PAMpal.accessors
#'
setMethod('[<-', 'AcousticStudy', function(x, i, value) {
    x@events[i] <- value
    validObject(x)
    x
})

#' @export
#' @rdname PAMpal.accessors
#'
setMethod('$', 'AcousticStudy', function(x, name) {
    '[['(x@events, name)
})

#' @export
#' @rdname PAMpal.accessors
#'
setMethod('$<-', 'AcousticStudy', function(x, name, value) {
    x@events[[name]] <- value
    validObject(x)
    x
})

#' @export
#' @rdname PAMpal.accessors
#'
setMethod('[[', signature(x='AcousticStudy', i='ANY', j='ANY'), function(x, i) {
# setMethod('[[', 'AcousticStudy', function(x, i, ...) {
    '[['(x@events, i)
})

#' @export
#' @rdname PAMpal.accessors
#'
setMethod('[[<-', 'AcousticStudy', function(x, i, value) {
    x@events[[i]] <- value
    validObject(x)
    x
})

#' @importFrom utils .DollarNames
#' @export
#'
.DollarNames.AcousticStudy <- function(x, pattern='') {
    grep(pattern, names(events(x)), value=TRUE)
}
