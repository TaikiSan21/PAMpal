## ---- PAMpalSettings Class ----------------------------------------------------
#' @title \code{PAMpalSettings} Class
#' @description An S4 class that stores settings related to all processing and analysis steps
#' done in PAMpal. A PAMpalSettings object will be the main input to any major function
#' in the PAMpal package.
#'
#' @slot db the full path to a PamGuard database file
#' @slot binaries a list with items "folder" containing the directory of the
#'   PamGuard binary files, and "list" containing the full path to each individual
#'   binary file.
#' @slot functions a named list of functions to apply to data read in by PAMpal.
#'   Should be named by the PamGuard module the function should be applied to.
#'   Currently supports "ClickDetector", "WhistlesMoans", and "Cepstrum".
#' @slot calibration a named list of calibration functions to apply while
#'   applying functions from the "functions" slot. Should named by the
#'   PamGuard module, same as the "functions"
#'
#' @author Taiki Sakai \email{taiki.sakai@@noaa.gov}
#' @export
#'
setClass('PAMpalSettings',
         slots = c(
             db = 'character',
             binaries = 'list',
             functions = 'list',
             calibration = 'list'
         ),
         prototype = prototype(
             db = character(0),
             binaries = list('folder'=character(0), 'list'=character(0)),
             functions = list('ClickDetector'=list(), 'WhistlesMoans'=list(), 'Cepstrum'=list()),
             calibration = list('ClickDetector'=list())
         )
)

setValidity('PAMpalSettings',
            function(object) {
                valid <- TRUE
                if(!all(c('folder', 'list') %in% names(object@binaries))) {
                    valid <- FALSE
                    cat('slot binaries must have items "folder" and "list"\n')
                }
                if(!all(c('ClickDetector', 'WhistlesMoans', 'Cepstrum') %in% names(object@functions))) {
                    valid <- FALSE
                    cat('slot functions must have items "ClickDetector", "WhistlesMoans", and "Cepstrum"\n')
                }
                valid
            }
)

#' @importFrom utils str
#'
setMethod('show', 'PAMpalSettings', function(object) {
    nBin <- length(object@binaries$list)
    nBinDir <- length(object@binaries$folder)
    nDb <- length(object@db)
    nCal <- length(object@calibration$ClickDetector)
    cat('PAMpalSettings object with:\n')
    cat(nDb, 'database(s)')
    if(nDb > 0) {
        # showDb <- basename(object@db)
        # if(nDb > 6) {
        #     showDb <- c(showDb[1:6], paste0('... (', nDb-6, ' more not shown)'))
        # }
        # cat(':\n ', paste(showDb, collapse='\n  '))
        cat(':\n ', printN(basename(object@db), n=6, collapse='\n  '))
    }
    cat('\n', nBinDir, ' binary folder(s) ', sep = '')
    if(nBinDir > 0) {
        cat('containing', nBin, 'binary files\n')
    } else {
        cat('\n')
    }
    # Print function names and args for each module
    for(m in seq_along(object@functions)) {
        cat(length(object@functions[[m]]), ' function(s) for module type "',
            names(object@functions)[m], '"\n', sep = '')
        for(f in seq_along(object@functions[[m]])) {
            cat(' "', names(object@functions[[m]])[f], '"\n  ', sep = '')
            cat(str(object@functions[[m]][[f]]))
        }
    }
    cat(nCal, 'click calibration function(s)\n')
})

#' Check if an Object is a PAMpalSettings
#'
#' Function to check if an object is a PAMpalSettings
#'
#' @param x object to check
#'
#' @export
#'
is.PAMpalSettings <- function(x) {
    inherits(x, 'PAMpalSettings')
}

## ---- AcousticEvent Class ---------------------------------------------------

# Acoustic event (obj) <--- this is really a list of AcEv? These need an ID for banter \acousticEvents
# Detector - named list [[detector name]] of lists    \\detector
# Data.table of detections w/ id
# possible image
# Localization - named list[[loc. type name]]         \\localization
# Data frame of positions
# Data Collection / Array Settings (obj)              \\settings
# Hydro sens, sample rate, whatever. Make an object and we figure out what it needs
# Visual data (obj)                                   \\visData
# Detection time, spp IDs, group size est, effort status. Multiple ways to read
# Behavioral (lul)                                    \\behavior
# erddap                                               \\erddap
# https://github.com/rmendels/Talks/blob/master/netCDF_Presentation/netcdf_opendap_erddap.Rmd
# Species classification - list of classifier objects \\species
# Method, prediction, assignment probabilities
# Duration? Files used? ID?

# setClassUnion('VisOrNULL', c('VisObsData', 'NULL'))

#' @title \code{AcousticEvent} Class
#' @description An S4 class storing acoustic detections from an Acoustic Event
#'   as well as other related metadata
#'
#' @slot id unique id or name for this event
#' @slot detectors a list of data frames that have acoustic detections and
#'   any measurements calculated on those detections. Each data frame is named
#'   by the detector that made the detection
#' @slot localizations a named list storing localizations, named by method
#' @slot settings a list of recorder settings
#' @slot species a list of species classifications for this event, named by
#'   classification method (ie. BANTER model, visual ID)
#' @slot files a list of files used to create this object, named by the type of
#'   file (ie. binaries, database)
#' @slot ancillary a list of miscellaneous extra stuff. Store whatever you want here
#'
#' @author Taiki Sakai \email{taiki.sakai@@noaa.gov}
#' @export
#'
setClass('AcousticEvent',
         slots = c(
             id = 'character',
             detectors = 'list',
             localizations = 'list',
             settings = 'list',
             species = 'list',
             files = 'list',
             ancillary = 'list'),
         prototype = prototype(id = character(), detectors=list(), localizations=list(),
                               settings=list(sr = NA_integer_, source = 'Not Found'), species=list(id=NA_character_),
                               files = list(), ancillary = list())
)

setValidity('AcousticEvent',
            function(object) {
                valid <- TRUE
                valid
            }
)
# Basic constructor
AcousticEvent <- function(id = character(), detectors=list(), localizations=list(),
                          settings=list(sr = NA_integer_, source = 'Not Found'),
                          species=list(id=NA_character_), files=list(), ancillary = list()) {
    new('AcousticEvent', id = id, detectors=detectors, localizations=localizations, settings=settings,
        species=species, files=files, ancillary = ancillary)
}

setMethod('show', 'AcousticEvent',
          function(object) {
              cat('AcousticEvent object "', id(object), '" with ',
                  length(object@detectors), ' detector(s): \n', sep='')
              cat(paste(names(object@detectors), collapse=', '))
          }
)

#' Check if an Object is an AcousticEvent
#'
#' Function to check if an object is an AcousticEvent
#'
#' @param x object to check
#'
#' @export
#'
is.AcousticEvent <- function(x) {
    inherits(x, 'AcousticEvent')
}

## ---- AcousticStudy Class ----------------------------------------------------------
# AcousticStudy class
# AcousticStudy (object)
# Files / folders (dbs, bins, vis, enviro)      \folders
# GPS                                           \gpsData
# Acoustic event (obj) <--- this is really a list of AcEv? These need an ID for banter \acousticEvents
# Detector - named list [[detector name]] of lists    \\detector
# Data.table of detections w/ id
# possible image
# Localization - named list[[loc. type name]]         \\localization
# Data frame of positions
# Data Collection / Array Settings (obj)              \\settings
# Hydro sens, sample rate, whatever. Make an object and we figure out what it needs
# Visual data (obj)                                   \\visData
# Detection time, spp IDs, group size est, effort status. Multiple ways to read
# Behavioral (lul)                                    \\behavior
# erddap                                               \\erddap
# Species classification - list of classifier objects \\species
# Method, prediction, assignment probabilities
# Detector settings - named list [[detector name]]   \detectorSettings
# Localization settings - named list [[ loc. type]]  \localizationSettings
# Some effort bullshit                               \effort
# ??????

# setOldClass(c('data.frame', 'data.table'))
#' @importClassesFrom data.table data.table
setClassUnion('dataframeORtable', members = c('data.frame', 'data.table'))

#' @title \code{AcousticStudy} Class
#' @description An S4 class storing acoustic data from an entire AcousticStudy
#' @slot id a unique id for the study
#' @slot events a list of \linkS4class{AcousticEvent} objects with
#'   detections from the AcousticStudy
#' @slot files a list of folders and files containing the AcousticStudy data
#' @slot gps a data frame of gps coordinates for the entire AcousticStudy
#' @slot pps the \linkS4class{PAMpalSettings} object used to create this object
#' @slot settings a named list of various settings for detectors, localizers, etc.
#' @slot effort something about effort lol
#' @slot models a place to store any models run on your data
#' @slot ancillary miscellaneous extra data
#'
#' @author Taiki Sakai \email{taiki.sakai@@noaa.gov}
#'
#' @importFrom utils packageVersion
#' @export
#'
setClass('AcousticStudy',
         slots = c(
             id = 'character',
             events = 'list',
             files = 'list',
             gps = 'dataframeORtable',
             pps = 'PAMpalSettings',
             settings = 'list',
             effort = 'dataframeORtable',
             models = 'list',
             ancillary = 'list'),
         prototype = prototype(
             id = character(),
             events=list(),
             files=list(db='None', binaries='None', visual='None', enviro='None'),
             gps=data.frame(),
             pps = new('PAMpalSettings'),
             settings=list(detectors=list(), localizations=list()),
             effort=data.frame(),
             models = list(),
             ancillary=list())
)

setValidity('AcousticStudy',
            function(object) {
                valid <- TRUE
                valid
            }
)

# Constructor
AcousticStudy <- function(id=NULL,
                          events=list(),
                          files=list(db=NA_character_, binaries=NA_character_, visual=NA_character_, enviro=NA_character_),
                          gps=data.frame(),
                          pps=new('PAMpalSettings'),
                          settings=list(detectors=list(), localizations=list()),
                          effort=data.frame(),
                          models = list(),
                          ancillary=list()) {
    if(is.null(id)) {
        id <- Sys.Date()
        message("No ID supplied for this AcousticStudy object, will use today's",
            ' date. Please assign a better name with id(study) <- "NAME"',
            '\nIn the future it is recommended to set the "id" argument.')
    }
    id <- as.character(id)
    fileTemp <- list(db=NA_character_, binaries=NA_character_, visual=NA_character_, enviro=NA_character_)
    for(n in names(files)) {
        fileTemp[[n]] <- files[[n]]
    }
    ancillary$version <- list(R = R.version.string,
                              PAMpal = packageVersion('PAMpal'))
    new('AcousticStudy', id=id, events=events, files=fileTemp, gps=gps,pps=pps,
        settings=settings, effort=effort, models=models, ancillary=ancillary)
}

#' Check if an Object is an AcousticStudy
#'
#' Function to check if an object is an AcousticStudy
#'
#' @param x object to check
#'
#' @export
#'
is.AcousticStudy <- function(x) {
    inherits(x, 'AcousticStudy')
}

setMethod('show', 'AcousticStudy',
          function(object) {
              cat('AcousticStudy object named ', id(object), ' with ',
                  length(events(object)), ' AcousticEvents.', sep='')
          }
)
