library(RNetCDF)
library(PAMpal)

renameSoundscopeRecordings <- function(nc, rec) {
    con <- open.nc(nc, write=TRUE)
    on.exit(close.nc(con))
    oldDirs <- var.get.nc(con, variable='audio_file_dir')
    if(all(dir.exists(unique(oldDirs)))) {
        cat('No change needed, all folders exist')
        return(invisible(oldDirs))
    }
    rec <- normalizePath(rec, winslash = '/')
    newDirs <- list.dirs(rec, full.names=TRUE, recursive=TRUE)
    
    newValues <- PAMpal:::fileMatcher(old=oldDirs, new=newDirs)
    changed <- newValues != oldDirs
    cat(sum(changed), 'out of', length(changed), 'values were changed.')
    var.put.nc(con, variable='audio_file_dir', data=newValues, start=1, count=length(newValues))
    invisible(newValues)
}
# change these paths so they are accurate for wherever youre working
baseDir <- 'PassiveAcoustics/DETECTOR_OUTPUT/PAMGUARD_ODONTOCETE_FIXED/Mid-Atlantic test/Spermwhale_Dolphin_Combo/'
recDir <- 'PassiveAcoustics_Soundfiles/BOTTOM_MOUNTED/NEFSC_VA/NEFSC_VA_202206_CB01/6102_48kHz/'

ncFile <- file.path(baseDir, 'testSoundscopeWhistles.nc')
# this will change paths in the netCDF to match above, should then be able to open
renameSoundscopeRecordings(ncFile, recDir)