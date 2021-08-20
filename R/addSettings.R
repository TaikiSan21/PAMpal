#' @title Add Settingss to a PAMpalSettings Object
#'
#' @description Adds settings to a PAMpalSettings object, usually
#'   from an XML file created by Pamguard's "Export XML Configuration"
#'
#' @param pps a \linkS4class{PAMpalSettings} object to add settings to
#' @param settings settings to add, either an XML file or a 
#' @param type one of 'xml' or 'list' indicating type of settings to add
#' @param verbose logical flag to show messages
#'
#' @return the same \linkS4class{PAMpalSettings} object as pps, with a new
#'   list of settings replacing what was previously in the "settings" slot
#'
#' @author Taiki Sakai \email{taiki.sakai@@noaa.gov}
#'
#' @examples
#'
#' # not recommended to create PPS like this, for example only
#' pps <- new('PAMpalSettings')
#' xmlSettings <- system.file('extdata', 'Example.xml', package='PAMpal')
#' pps <- addSettings(pps, xmlSettings, type='xml')
#'
#' @importFrom tcltk tk_choose.files
#' @export
#'
addSettings <- function(pps, settings=NULL, type=c('xml', 'list'), verbose=TRUE) {
  type <- match.arg(type)
  if(is.PAMpalSettings(settings)) {
    settings <- settings@settings
    type <- 'list'
  }
  if(is.null(settings)) {
    cat('Please select the XML settings file.\n')

    settings <- tk_choose.files(caption = 'Choose an XML Settings File:',
                              default = getwd(), multi=FALSE)
    type <- 'xml'
    if(is.null(settings) || is.na(settings)) {
      message('No file chosen')
      return(pps)
    }
  }
  # Case when cancelled, dont error
  switch(type,
         'list' = {
           if(validSettings(settings)) {
             pps@settings <- settings
           }
         },
         'xml' = {
           if(!file.exists(settings)) {
             warning(paste0('XML settings file ', settings, ' does not exist'))
             return(pps)
           }
           settings <- loadPamguardXML(settings)
           pps@settings <- settings
         })
  if(verbose) {
    cat(paste0('Adding settings for ', length(settings$detectors), ' detectors.\n'))
  }
  pps
}

validSettings <- function(x) {
  is.list(x) && all(c('sources', 'detectors') %in% names(x))
}