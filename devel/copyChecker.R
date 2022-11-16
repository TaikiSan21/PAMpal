# original: base folder of original data
# copy: base folder of copied data
# progress: show progress bar or not
# usage: copyChecker('./Data/Original', './Data/Copy')
# may take a long time or crash with huge folders
copyChecker <- function(original, copy, progress=TRUE) {
    oldFiles <- list.files(original, full.names = FALSE, recursive = TRUE)
    newFiles <- list.files(copy, full.names = FALSE, recursive = TRUE)
    missing <- !(oldFiles %in% newFiles)
    if(progress) {
        pb <- txtProgressBar(min=0, max=length(oldFiles[!missing])+length(newFiles), style=3)
        ix <- 1
    }
    oldSize <- unlist(sapply(oldFiles[!missing], function(x) {
        if(progress) {
            setTxtProgressBar(pb, value=ix)
            ix <<- ix +1
        }
        file.size(file.path(original, x))
    }))
    names(oldSize) <- NULL
    newSize <- unlist(sapply(newFiles, function(x) {
        if(progress) {
            setTxtProgressBar(pb, value=ix)
            ix <<- ix +1
        }
        file.size(file.path(copy, x))
    }))
    names(newSize) <- NULL
    sizeDiff <- oldSize - newSize
    changed <- sizeDiff != 0
    cat('\n', sum(missing), ' missing files\n',
        sum(abs(sizeDiff)), ' size difference (bytes) between matching files',
        sep='')
    list(missing=oldFiles[missing],
         changed=oldFiles[changed],
         oldSize=oldSize,
         newSize=newSize)
}
