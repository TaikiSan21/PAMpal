#' @title Add a Function to a PAMpalSettings Object
#'
#' @description Adds a new function to the "function" slot in a PAMpalSettings
#'   object. Must be run interactively, user will be prompted to assign
#'   values for any parameters in the function to be added
#'
#' @param pps a \linkS4class{PAMpalSettings} object to add a function to
#' @param fun function to add OR another \linkS4class{PAMpalSettings} object.
#'   In this case all functions from the second object will be added to \code{pps}
#' @param module Pamguard module output this function should act on, one of
#'   ClickDetector, WhistlesMoans, Cepstrum, or GPLDetector. If \code{NULL} (default), user
#'   will be prompted to select which module it applies to
#' @param verbose logical flag to show messages
#' @param default logical flag to use default function parameters if present
#' @param \dots named arguments to pass to function being added
#'
#' @return the same \linkS4class{PAMpalSettings} object as pps, with the function
#'   \code{fun} added to the "functions" slot
#'
#' @author Taiki Sakai \email{taiki.sakai@@noaa.gov}
#'
#' @examples
#'
#' # not recommended to create a pps like this, for example only
#' pps <- new('PAMpalSettings')
#' if(interactive()) pps <- addFunction(pps, standardClickCalcs)
#' pps <- addFunction(pps, roccaWhistleCalcs, module='WhistlesMoans')
#'
#' @importFrom utils menu
#' @export
#'
addFunction <- function(pps, fun, module=NULL, verbose = TRUE, default=FALSE,  ...) {
    modsAllowed <- c('ClickDetector', 'WhistlesMoans', 'Cepstrum', 'GPLDetector')
    if(is.PAMpalSettings(fun)) {
        for(m in modsAllowed) {
            newFuns <- fun@functions[[m]]
            if(verbose) {
                cat('Adding', length(newFuns), 'functions to module type', m, '\n')
            }
            pps@functions[m] <- list(c(pps@functions[[m]], newFuns))
        }
        return(pps)
    }
    if(is.null(attr(fun, 'fname'))) {
        fname <- deparse(substitute(fun))
    } else {
        fname <- attr(fun, 'fname')
    }
    if(is.null(module) ||
       !(module %in% modsAllowed)) {
        chooseMod <- menu(choices = modsAllowed,
                          title = c('What module type is this function for?'))
        if(chooseMod==0) stop('You must set a module type for this function.')
        module <- modsAllowed[chooseMod]
    }
    if(verbose) {
        cat('Adding function "', fname, '":\n', sep = '')
    }
    oldnames <- names(pps@functions[[module]])
    fun <- functionParser(fun, default=default, ...)
    # function checker
    if(functionChecker(fun, module)) {
        pps@functions[module] <- list(c(pps@functions[[module]], fun))
        names(pps@functions[[module]]) <- c(oldnames, fname)
    } else {
        warning('Unable to add function ', fname, ', it did not create the expected output.')
    }
    pps
}

# I put a function in yo function cuz i heard you like functions
functionParser <- function(fun, skipArgs = c('data', 'calibration', '...'), default=FALSE, ...) {
    argList <- formals(fun)
    dotList <- list(...)
    dotArgs <- names(argList)[names(argList) %in% names(dotList)]
    for(d in dotArgs) {
        argList[d] <- dotList[d]
        skipArgs <- c(skipArgs, d)
    }
    if(isTRUE(default)) {
        for(arg in names(argList)) {
            # "name" class means no default value 
            if(inherits(argList[[arg]], 'name')) {
                next
            }
            skipArgs <- c(skipArgs, arg)
        }
    }
    toSet <- names(argList)[!(names(argList) %in% skipArgs)]
    if(length(toSet) > 0) {
        for(a in toSet) {
            cat('Set a value for parameter "', a, '", please put quotes around strings', sep='')
            if(inherits(argList[[a]], 'name')) {
                cat(' (no default value found):')
            } else if(inherits(argList[[a]], 'NULL')) {
                cat(' (default value is NULL):')
            } else if(inherits(argList[[a]], 'call')) {
                cat(' (default value is ', deparse(argList[[a]]), '):', sep='')
            } else {
                cat(' (default value is ', argList[[a]], '):', sep = '')
            }
            val <- readline()
            # If it evals properly, do that. Otherwise its prob a string so leave it
            newVal <- tryCatch({
                # empty so you dont accidentally grab vars, or do you want this????
                # just change to globalenv() if you want that
                eval(parse(text=val), envir = globalenv())
            },
            error = function(e) {
                val
            })
            if(is.null(newVal)) next
            argList[a] <- newVal
        }
    }
    formals(fun) <- argList
    fun
}

functionChecker <- function(fun, module) {
    switch(module,
           'ClickDetector' = clickChecker(fun),
           'WhistlesMoans' = whistleChecker(fun),
           'Cepstrum' = cepstrumChecker(fun),
           'GPLDetector' = gplChecker(fun),
           FALSE
    )
}

clickChecker <- function(fun) {
    good <- TRUE

    testThisClick <- try(fun(data=PAMpal::testClick))

    if(inherits(testThisClick, 'try-error')) {
        message('Click function did not run succesfully.')
        message('Error: ', attr(testThisClick, 'condition')$message)
        return(FALSE)
    }
    if(is.null(testThisClick)) {
        message('Click function returned nothing.')
        return(FALSE)
    }
    if(nrow(testThisClick) != 2) {
        message('Click functions should return 1 row for each channel.')
        good <- FALSE
    }
    good
}

whistleChecker <- function(fun) {
    good <- TRUE

    testThisWhistle <- try(fun(data=PAMpal::testWhistle))

    if(inherits(testThisWhistle, 'try-error')) {
        message('Whistle function did not run successfully.')
        message('Error: ', attr(testThisWhistle, 'condition')$message)
        return(FALSE)
    }
    if(is.null(testThisWhistle)) {
        message('Whistle function returned nothing.')
        return(FALSE)
    }
    if(nrow(data.frame(testThisWhistle)) != 1) {
        message('Whistle functions should return a single row for each contour.')
        return(FALSE)
    }
    good
}

cepstrumChecker <- function(fun) {
    good <- TRUE
    testThisCeps <- try(fun(data=PAMpal::testCeps))

    if(inherits(testThisCeps, 'try-error')) {
        message('Cepstrum function did not run successfully.')
        message('Error: ', attr(testThisCeps, 'condition')$message)
        return(FALSE)
    }
    if(is.null(testThisCeps)) {
        message('Cepstrum function returned nothing.')
        return(FALSE)
    }
    if(nrow(data.frame(testThisCeps)) != 1) {
        message('Cepstrum function should return a single row for each detection.')
        return(FALSE)
    }
    good
}

gplChecker <- function(fun) {
    good <- TRUE
    testThisGpl <- try(fun(data=PAMpal::testGPL))

    if(inherits(testThisGpl, 'try-error')) {
        message('GPL function did not run successfully.')
        message('Error: ', attr(testThisGpl, 'condition')$message)
        return(FALSE)
    }
    if(is.null(testThisGpl)) {
        message('GPL function returned nothing.')
        return(FALSE)
    }
    if(nrow(data.frame(testThisGpl)) != 1) {
        message('GPL functions should return a single row for each contour.')
        return(FALSE)
    }
    good
}