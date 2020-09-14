#' @title Remove Binaries from a PAMpalSettings Object
#'
#' @description Remove a binary folder and associated files from the "binaries"
#'   slot in a PAMpalSettings object.
#'
#' @param pps a \linkS4class{PAMpalSettings} object to remove binaries from
#' @param index index indicating which binary folders to remove. Can be a vector
#'   if you want to remove multiple folders. If missing user is prompted to
#'   select a folder from a list, will only show up to the first 20. You can
#'   easily remove all of the folders with a large index like \code{1:1000}
#'
#' @return the same \linkS4class{PAMpalSettings} object as pps, with the binary
#'   folders and files associated with those folders removed from the "binaries"
#'   slot.
#'
#' @author Taiki Sakai \email{taiki.sakai@@noaa.gov}
#'
#' @importFrom utils menu
#' @export
#'
removeBinaries <- function(pps, index=NULL) {
    if(is.null(index)) {
        if(length(pps@binaries$folder) > 20) {
            warning('Only showing first 20 binary folders.')
            choices <- pps@binaries$folder[1:20]
        }
        choices <- pps@binaries$folder
        index <- menu(title = 'Choose a folder to remove:',
                      choices = choices)
        if(index==0) return(pps)
    }
    if(max(index) > length(pps@binaries$folder)) warning('Index too large, no folder to remove.')
    dropNames <- pps@binaries$folder[index]
    for(f in dropNames) {
        if(is.na(f)) next
        pps@binaries$list <- pps@binaries$list[!grepl(f, pps@binaries$list, fixed=TRUE)]
    }
    pps@binaries$folder <- pps@binaries$folder[-index]
    pps
}
