#' @title Add a Calibration File to a PAMpalSettings Object
#'
#' @description Adds a new calibration function to the "calibration" slot
#'   of a PAMpalSettings object. Interactively asks user to supply
#'   file and other parameters if not supplied.
#'
#' @param pps a \linkS4class{PAMpalSettings} object to add a database to
#' @param calFile a calibration file name. Must be csv format with two columns.
#'   The first column must be the frequency (in Hz), and the second column must be the
#'   sensitivity (in dB), and the columns should be labeled \code{Frequency}
#'   and \code{Sensitivity}. Can also be supplied as a dataframe in which
#'   case the \code{calName} argument should also be set
#' @param module the Pamguard module type this calibration should be applied to,
#'   for now this is only for ClickDetector modules. This is left as an option
#'   for future-proofing purposes but should not be needed.
#' @param calName the name to assign to the calibration function, defaults to
#'   the file name and only needs to be set if supplying a dataframe instead of
#'   a csv file
#' @param all logical flag whether or not to apply calibration to all functions
#'   without asking individually, recommended to stay as \code{FALSE}
#' @param units a number from 1 to 3 specifying the units of the calibration file,
#'   number corresponds to dB re V/uPa, uPa/Counts, or uPa/FullScale respectively.
#'   A NULL (default) or other value will prompt user to select units.
#'
#' @return the same \linkS4class{PAMpalSettings} object as pps, with the calibration
#'   function added to the \code{calibration} slot.
#'
#' @details When adding a calibration, you will be asked what units your calibration
#'   value is in. The wave clips stored by Pamguard are values from -1 to 1, so if
#'   your calibration is expecting different units then this needs to be accounted
#'   for in order to get an accurate SPL value. For V / uPa you must know the voltage
#'   range of your recording equipment, and for calibrations expecting Count data you
#'   must know the bit rate of your recordings. If your calibration is already relative
#'   to full-scale, then nothing needs to be adjusted. If you don't know the units of
#'   your calibration and you are only interested in relative dB levels, then you can
#'   select the Full-Scale options.
#'
#'   The calibration function created takes frequency (in Hz) as input and
#'   outputs the associated dB value that needs to be added to correct the
#'   power spectrum of a signal. If the input is a matrix or dataframe, the first
#'   column is assumed to be frequency.
#'
#' @author Taiki Sakai \email{taiki.sakai@@noaa.gov}
#'
#' @examples
#'
#' pps <- new('PAMpalSettings')
#' calFile <- system.file('extdata', 'calibration.csv', package='PAMpal')
#' pps <- addCalibration(pps, calFile, all = TRUE, units=3)
#' calClick <- function(data, calibration=NULL) {
#'     standardClickCalcs(data, calibration=calibration, filterfrom_khz = 0)
#' }
#' pps <- addFunction(pps, calClick, module = 'ClickDetector')
#' pps <- applyCalibration(pps, all=TRUE)
#' pps
#'
#' @importFrom utils read.csv menu
#' @importFrom gam gam s predict.Gam
#' @importFrom seewave spec
#' @importFrom tcltk tk_choose.files
#' @export
#'
addCalibration <- function(pps, calFile=NULL, module='ClickDetector', calName=NULL, all=FALSE, units=NULL) {
    if(is.PAMpalSettings(calFile)) {
        funsToAdd <- calFile@calibration[[module]]
        for(i in seq_along(funsToAdd)) {
            pps <- addCalibration(pps,
                                  get('CAL', envir = environment(calFile@calibration[[module]][[i]])),
                                  module=module,
                                  calName=names(calFile@calibration[[module]])[i])
        }
        return(pps)
    }
    modsAllowed <- c('ClickDetector')
    if(!(module %in% modsAllowed)) {
        warning("I don't know how to add calibration for module type", module)
        return(pps)
    }
    if(is.null(calFile)) {
        cat('Please choose a calibration file (must be csv format).\n')
        calFile <- tk_choose.files(caption = 'Select calibration file:')
    }
    if(length(calFile) == 0) return(pps)

    calFun <- makeCalibration(calFile, units)
    if(is.null(calFun)) return(pps)

    oldNames <- names(pps@calibration[[module]])
    if(is.null(calName)) {
        calName <- ifelse(is.character(calFile), calFile, 'Calibration')
    }

    pps@calibration[module] <- list(c(pps@calibration[[module]], calFun))

    names(pps@calibration[[module]]) <- c(oldNames, calName)
    pps <- applyCalibration(pps, module=module, all=all)
    pps
}

# Do checks, return calibration as a function
makeCalibration <- function(calFile, units=NULL) {

    # POSSIBLY add at top single value calibration

    if(is.character(calFile)) {
        if(!grepl('csv', calFile)) {
            warning('Calibration file must be provided in csv format with columns',
                    '"Frequency" and "Sensitivity", calibration has not been added.')
            return(NULL)
        }
        CAL <- read.csv(calFile, header = FALSE, stringsAsFactors = FALSE)
        isRowname <- is.na(CAL[1, ])
        if(any(isRowname)) {
            warning('First column appears to be row names or ids, these are being removed.')
            CAL <- CAL[, !isRowname]
        }
        hasHeader <- grepl('[A-z]', CAL[1, 1]) && grepl('[A-z]', CAL[1, 2])
        if(hasHeader) {
            CAL <- CAL[-1, ]
        }
        CAL[, 1] <- suppressWarnings(as.numeric(CAL[, 1]))
        CAL[, 2] <- suppressWarnings(as.numeric(CAL[, 2]))
        if(any(is.na(CAL[, 1])) ||
           any(is.na(CAL[, 2]))) {
            warning('Non-numeric values found in calibration file, please check.',
                    'Calibration has not been applied.')
            return(NULL)
        }
    }
    if(is.data.frame(calFile)) {
        CAL <- calFile

    }
    colnames(CAL) <- c('Frequency', 'Sensitivity')
    CAL$Frequency <- as.numeric(CAL$Frequency)
    if(max(CAL$Frequency) < 1e3) {
        warning('It appears that frequency was supplied in kHz, converting to Hz')
        CAL$Frequency <- CAL$Frequency * 1e3
    }
    CAL$Sensitivity <- as.numeric(CAL$Sensitivity)
    # Convert to always dB in numerator, this should always work !!NVM doing this after predict
    # should be better so plotting shows negative like they supplied so no confusing?
    # if(median(CAL$Sensitivity, na.rm = TRUE) < 0) {
    #     CAL$Sensitivity <- CAL$Sensitivity * -1
    # }
    unitOpts <- paste0('dB re ', c('V/uPa', 'uPa/Counts', 'uPa/FullScale'))
    if(is.null(units) ||
       !(units %in% 1:3)) {
        unitChoice <- menu(title = paste0('What are the units of your calibration?\n',
                                      "(If you only need relative dB values choose the FullScale option)"),
                       choices = unitOpts)
    } else {
        unitChoice <- units
    }
    if(unitChoice == 0) {
        stop('No unit selection made, stopping calibration.')
    }
    if(unitChoice == 1) {
        scaleRead <- readline(prompt='What is your voltage range?')
        SCALE <- suppressWarnings(as.numeric(scaleRead))
    }
    if(unitChoice == 2) {
        scaleRead <- readline(prompt='What is the bit rate of your data?')
        SCALE <- suppressWarnings(2 ^ (as.numeric(scaleRead) - 1))
    }
    if(unitChoice == 3) {
        SCALE <- 1
    }
    if(is.na(SCALE)) {
        stop('Unable to convert input to a numeric value, stopping calibration.')
    }

    # need to convert all to assumed dB / FullScale, which is what corresponds to converting
    # -1 to 1 range to actual dB

    # say that fullscale or end-end option is if you dont care about absolute only relative
    # other require input numbers

    # If Volts multiple spec by V, if counts multiple spec by 2^Bit
    CALMODEL <- suppressWarnings(gam(Sensitivity ~ s(Frequency, 50), data = CAL))

    # Applies to power spectrum 20* log10(spec).
    # Usually subtracted? Added?
    # ORIGINAL EXAMPLE FROM JAY HAS NEGATIVE VALUES IN CAL FILE ~-170

    # spec has column 1 frequency kHz b/c thats what seewave::spec outputs, column 2 dB in 20*log10
    calFun <- function(freq, ...) {
        # if take in spec we add 20*log10(SCALE)
        # specData <- spec(wave = wave * SCALE, f = sr, norm = FALSE,
        #                           correction = 'amplitude', plot = FALSE, ...)
        # # This is for integer overflow spec jankiness fixing
        # freq <- seq(from = 0, by = specData[2, 1] - specData[1,1], length.out = nrow(specData))
        # specData <- data.frame(Frequency = freq * 1e3,
        #                        Sensitivity = specData[, 2])
        # specData$Sensitivity <- 20 * log10(specData$Sensitivity)
        # specData$Sensitivity[!is.finite(specData$Sensitivity)] <- NA
        # specData <- data.frame(Frequency = spec[, 1] * 1e3,
        #                        Sensitivity = spec[, 2])
        # specData$Sensitivity  <- specData$Sensitivity + 20*log10(SCALE)
        if(is.matrix(freq) ||
           is.data.frame(freq)) {
            freq <- freq[, 1]
        }
        predicted <- predict.Gam(CALMODEL, newdata = data.frame(Frequency = freq))
        if(median(predicted) < 0) {
            predicted <- predicted * -1
        }
        # predicted <- predicted - predicted[1]
        # this is + or - WHO KNOWS WTF
        # specData$Sensitivity <- specData$Sensitivity + predicted
        # specData$Sensitivity <- specData$Sensitivity - max(specData$Sensitivity)
        # colnames(specData) <- c('Frequency', 'dB')
        # possibly i should save function and values?
        predicted + 20*log10(SCALE)
    }
    class(calFun) <- c('calibration', 'function')
    calFun
}

#' @export
#'
plot.calibration <- function(x, ...) {
    mdl <- get('CALMODEL', envir = environment(x))
    plot(mdl)
}

# Add calibration name to click funs
#' @rdname addCalibration
#' @export
#'
applyCalibration <- function(pps, module = 'ClickDetector', all=FALSE) {
    calList <- pps@calibration[[module]]
    if(length(calList) == 0) {
        message('No calibration functions found in PPS.\n')
        return(pps)
    }
    if(length(calList) == 1) {
        calName <- names(calList)
        calFun <- calList[[1]]
    }
    if(length(calList) > 1) {
        chooseFun <- menu(choices = names(calList),
                          title = 'Which calibration function should be applied?')
        if(chooseFun == 0) {
            warning('You must choose a function.')
            return(pps)
        }
        calName <- names(calList)[chooseFun]
        calFun <- calList[[chooseFun]]
    }
    argList <- lapply(pps@functions[[module]], formals)
    hasCal <- sapply(argList, function(x) 'calibration' %in% names(x))
    if(length(hasCal) == 0) {
        message('Could not find any functions with "calibration" as an input.\n')
        return(pps)
    }
    hasCal <- which(hasCal)
    for(i in hasCal) {
        if(!all) {
            yn <- menu(choices = c('Yes', 'No'),
                       title = paste0('Should we apply calibration ', calName,
                                      ' to function ', names(argList)[i], '?'))
            if(yn == 2 || yn == 0) next
        }
        formals(pps@functions[[module]][[i]])['calibration'] <- calName
    }
    pps
}

# do janky environment search to actually get functiont
findCalibration <- function(calName, module = 'ClickDetector') {
    allGlobal <- ls(envir = .GlobalEnv)
    allPps <- sapply(allGlobal, function(x) {
        obj <- get(x, envir = .GlobalEnv)
        if('PAMpalSettings' %in% class(obj)) {
            obj
        } else if(is.AcousticStudy(obj)) {
            obj@pps
        } else {
            NULL
        }
        })
    allPps <- allPps[!sapply(allPps, is.null)]
    hasCal <- which(sapply(allPps, function(x) calName %in% names(x@calibration[[module]])))
    if(length(hasCal) == 0) {
        stop('Could not find calibration function named ', calName,
             ' please re-apply calibration.')
    }
    if(length(hasCal) > 1) {
        # warning('Found more than one instance of calibration function named ', calName,
        #         ', using the first one.')
        useCal <- hasCal[1]
    } else {
        useCal <- hasCal[1]
    }
    allPps[[useCal]]@calibration[[module]][[calName]]
}
