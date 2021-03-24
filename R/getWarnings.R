#' Get Warning Messages
#' 
#' Accessor to easily get all warning messages for \code{x}
#' 
#' @param x an \linkS4class{AcousticStudy} or \linkS4class{AcousticEvent}
#'   object
#' 
#' @return a list of warning messages, named by the function call that created
#'   the warning
#'
#' @author Taiki Sakai \email{taiki.sakai@@noaa.gov}
#' 
#' @export
#' 
getWarnings <- function(x) {
  ancillary(x)$warnings
}

#' @importFrom PAMmisc squishList
#' 
pamWarning <- function(..., which=1, n=6) {
  warnState <- options('pamWarn')[['pamWarn']]
  if(!is.null(warnState) &&
     warnState == -1) {
    return(TRUE)
  }
  dotList <- list(...)
  showMessage <- paste0(
    sapply(dotList, function(x) printN(x, n)), collapse=''
  )
  warning(showMessage, call.=FALSE)
  saveMessage <- paste0(
    sapply(dotList, function(x) printN(x, Inf)), collapse=''
  )
  oldWarn <- mget('PAMWARNING', envir = sys.frame(1), ifnotfound=NA)[['PAMWARNING']]
  if(which > sys.nframe()) {
    which <- sys.nframe() - 1
  }
  caller <- deparse(sys.call(which)[[1]])
  newWarn <- data.frame(time=Sys.time(), functionName=caller, message=saveMessage)
  
  if(is.data.frame(oldWarn)) {
    newWarn <- rbind(oldWarn, newWarn)
  } 
  assign('PAMWARNING', value = newWarn, envir = sys.frame(1))
  invisible(TRUE)
}

.addPamWarning <- function(x) {
  warns <- mget('PAMWARNING', envir = sys.frame(1), ifnotfound = NA)[['PAMWARNING']]
  if(length(warns) == 1 &&
     is.na(warns)) {
    # warns <- data.frame(time=Sys.time(), functionName=deparse(sys.call(1)[[1]]), message='No warnings')
    return(x)
  }
  # for(i in seq_along(warns)) {
  #   warns[[i]] <- unique(warns[[i]])
  # }
  ancillary(x)$warnings <- tibble(rbind(ancillary(x)$warnings, warns))
  x
}

suppressPamWarnings <- function(expr) {
  oldState <- options('pamWarn')[['pamWarn']]
  options(pamWarn=-1)
  on.exit(options(pamWarn=oldState))
  suppressWarnings(expr)
}
