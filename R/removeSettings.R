#' @title Remove Settings from a PAMpalSettings Object
#'
#' @description Remove all settings from the "settings" slot in a PAMpalSettings
#'   object.
#'
#' @param pps a \linkS4class{PAMpalSettings} object to remove settings from
#' 
#' @return the same \linkS4class{PAMpalSettings} object as pps, with all
#'   settings removed from the "settings" slot
#'
#' @author Taiki Sakai \email{taiki.sakai@@noaa.gov}
#'
#' @examples
#'
#' exPps <- new('PAMpalSettings')
#' exPps <- addSettings(exPps, system.file('extdata', 'Example.xml', package='PAMpal'))
#' removeSettings(exPps)
#'
#' @export
#'
removeSettings <- function(pps) {
  pps@settings <- list()
  pps
}
