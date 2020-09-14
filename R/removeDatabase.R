#' @title Remove a Database from a PAMpalSettings Object
#'
#' @description Remove a database from the "db" slot in a PAMpalSettings
#'   object.
#'
#' @param pps a \linkS4class{PAMpalSettings} object to remove a database from
#' @param index index indicating which database(s) to remove. Can be a vector
#'   if you want to remove multiple databases. If missing user is prompted to
#'   select a database from a list, will only show up to the first 20. You can
#'   easily remove all of the databases with a large index like \code{1:1000}
#'
#' @return the same \linkS4class{PAMpalSettings} object as pps, with the database(s)
#'   removed from the "db" slot
#'
#' @author Taiki Sakai \email{taiki.sakai@@noaa.gov}
#'
#' @importFrom utils menu
#' @export
#'
removeDatabase <- function(pps, index=NULL) {
    if(is.null(index)) {
        choices <- pps@db
        if(length(choices) > 20) {
            warning('Only showing first 20 databases.')
            choices <- choices[1:20]
        }
        index <- menu(title = 'Choose a database to remove:',
                      choices = choices)
        if(index==0) return(pps)
    }
    if(max(index) > length(pps@db)) warning('Index too large, no database to remove.')
    pps@db <- pps@db[-index]
    pps
}
