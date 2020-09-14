#' @title Add Binaries to a PAMpalSettings Object
#'
#' @description Adds a new function to the "function" slot in a PAMpalSettings
#'   object.
#'
#' @param pps a \linkS4class{PAMpalSettings} object to add a database to
#' @param binFolder a folder of binaries to add
#'
#' @return the same \linkS4class{PAMpalSettings} object as pps, with the binary
#'   files contained in \code{binFolder} added to the "binaries" slot. Only
#'   binary files for Click Detector and WhistlesMoans modules will be added,
#'   since these are the only types PAMpal currently knows how to process
#'   (last updated v 0.7.0)
#'
#' @author Taiki Sakai \email{taiki.sakai@@noaa.gov}
#'
#' @examples
#'
#' # not recommended to create PPS like this, for example only
#' pps <- new('PAMpalSettings')
#' binFolder <- system.file('extdata', 'Binaries', package='PAMpal')
#' pps <- addBinaries(pps, binFolder)
#' pps
#'
#' @importFrom tcltk tk_choose.dir
#' @export
#'
addBinaries <- function(pps, binFolder=NULL) {
    binList <- NULL
    if(is.PAMpalSettings(binFolder)) {
        binList <- binFolder@binaries$list
        binFolder <- binFolder@binaries$folder
        exists <- file.exists(binList)
        if(any(!exists)) {
            binList <- binList[exists]
            cat(sum(!exists), 'binary files did not exist.\n')
        }
    }
    if(is.null(binFolder)) {
        cat('Please select the folder where the binaries are stored.\n')
        # if(!rstudioapi::isAvailable('1.1.287')) {
        #     binFolder <- choose.dir(caption = 'Choose Binary Folder:')
        # } else {
        #     binFolder <- rstudioapi::selectDirectory(caption = 'Choose Binary Folder:', path = getwd())
        # }
        binFolder <- tk_choose.dir(caption = 'Choose Binary Folder:',default = getwd())
    }
    # Case when cancelled, dont error
    if(is.null(binFolder) || is.na(binFolder)) {
        cat('No folder chosen')
        return(pps)
    }
    if(!dir.exists(binFolder)) {
        cat(paste0('Binary folder ', binFolder, ' does not exist'))
        return(pps)
    }
    pps@binaries$folder <- unique(c(pps@binaries$folder, binFolder))
    if(is.null(binList) ||
       length(binList) == 0) {
        cat('Getting list of all binary files in folder. This may take a while...\n')
        # only have functions for Clicks & Whistles right now, filter out so we dont get garbage
        # warning overflow later
        binList <- list.files(binFolder, recursive = TRUE, full.names = TRUE, pattern ='(Clicks|WhistlesMoans).*pgdf$')
    }
    cat('Adding', length(binList), 'binary files from', length(binFolder), 'folders\n')
    pps@binaries$list <- unique(c(pps@binaries$list, binList))
    pps
}
