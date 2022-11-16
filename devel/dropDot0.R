dropDot0 <- function(dir=NULL) {
    if(is.null(dir)) {
        dir <- choose.dir()
    }
    files <- list.files(dir, pattern='\\.0\\.wav$', full.names=TRUE)
    if(length(files) == 0) {
        cat('No files found to change.')
        return(invisible(TRUE))
    }
    changed <- sapply(files, function(x) {
        if(!file.exists(x)) {
            message('File ', x, ' does not exist')
            return(FALSE)
        }
        canWrite <- file.access(x, 2) == 0
        if(!canWrite) {
            message('No write access for file ', x)
            return(FALSE)
        }
        newName <- gsub('\\.0\\.wav$', '.wav', x)
        message('Trying to change file ', basename(x), ' to ',
                basename(newName))
        file.rename(from=x, to=newName)
    })
    cat('Changed ', sum(changed), ' out of ', length(changed), ' files.', sep='')
    invisible(changed)
}