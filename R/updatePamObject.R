#' Update PAMpal S4 Object
#' 
#' Updates older versions of PAMpal's S4 objects to stop "validOjbect"
#' warning messages
#' 
#' @details As of v0.12.0 this updates any previous version's PAMpalSettings
#'   objects to have the new "settings" slot, as well as updating any 
#'   PAMpalSettings objects within an AcousticStudy
#'   
#' @param x an \linkS4class{AcousticStudy}, \linkS4class{AcousticEvent},
#'   or \linkS4class{PAMpalSettings} object
#' 
#' @return the same object as \code{x} with any slot changes made
#' 
#' @author Taiki Sakai \email{taiki.sakai@@noaa.gov}
#' 
#' @examples
#' 
#' \dontrun{
#' pps <- new('PAMpalSettings')
#' # manually breaking this S4 class, don't try this at home
#' attr(pps, 'settings') <- NULL
#' # This will now give an error
#' pps
#' pps <- updatePamObject(pps)
#' # Fixed!
#' pps
#' }
#' 
#' @export
#' 
updatePamObject <- function(x) {
  if(is.AcousticEvent(x)) {
    return(x)
  }
  if(is.AcousticStudy(x)) {
    x@pps <- updatePamObject(x@pps)
    return(x)
  }
  if(is.PAMpalSettings(x)) {
    if(is.null(attr(x, 'settings'))) {
      x@settings <- list(file=NA_character_,
                         sources=list(),
                         detectors=list(),
                         raw=NA_character_)
    }
    return(x)
  }
}
