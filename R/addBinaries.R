#' @title Add Binaries to a PAMpalSettings Object
#'
#' @description Adds more binary files to the "binaries" slot of a
#'   PAMpalSettings object. Interactively asks user to supply
#'   folder location if not provided.
#'
#' @param pps a \linkS4class{PAMpalSettings} object to add binary files to
#' @param folder a folder of binaries to add, all subfolders will also be added
#' @param verbose logical flag to show messages
#'
#' @return the same \linkS4class{PAMpalSettings} object as pps, with the binary
#'   files contained in \code{folder} added to the "binaries" slot. Only
#'   binary files for Click Detector and WhistlesMoans modules will be added,
#'   since these are the only types PAMpal currently knows how to process
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
addBinaries <- function(pps, folder=NULL, verbose=TRUE) {
    binList <- NULL
    if(is.PAMpalSettings(folder)) {
        binList <- folder@binaries$list
        folder <- folder@binaries$folder
        exists <- file.exists(binList)
        if(any(!exists)) {
            binList <- binList[exists]
            if(verbose) {
                cat(sum(!exists), 'binary files did not exist.\n')
            }
        }
    }
    if(is.null(folder)) {
        cat('Please select the folder where the binaries are stored.\n')
        folder <- tk_choose.dir(caption = 'Choose Binary Folder:',default = getwd())
    }
    # Case when cancelled, dont error
    if(all(is.null(folder)) || all(is.na(folder))) {
        message('No folder chosen')
        return(pps)
    }
    foldExist <- sapply(folder, dir.exists)
    if(!all(foldExist)) {
        warning(paste0('Binary folder ', printN(folder[!foldExist]), ' does not exist'))
        folder <- folder[foldExist]
    }
    if(length(folder) == 0) {
        message('No folders to add binaries from.')
        return(pps)
    }
    folder <- normalizePath(folder, winslash = '/')
    pps@binaries$folder <- unique(c(pps@binaries$folder, folder))
    if(is.null(binList) ||
       length(binList) == 0) {
        if(verbose) {
            cat('Getting list of all binary files in folder. This may take a while...\n')
        }
        # only have functions for Clicks & Whistles right now, filter out so we dont get garbage
        # warning overflow later
        for(f in folder) {
            binList <- c(binList,
                         list.files(f, recursive = TRUE, full.names = TRUE, pattern = ppVars()$binPattern))
        }
    }
    binList <- binList[!grepl('GPL_State', binList)]
    if(verbose) {
        cat('Adding', length(binList), 'binary files from', length(folder), 'folders\n')
    }
    pps@binaries$list <- unique(c(pps@binaries$list, binList))
    pps
}
