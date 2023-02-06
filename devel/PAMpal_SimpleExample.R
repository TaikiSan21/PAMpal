# PAMpal simple example 
# Its on CRAN. Yay!
# install.packages('PAMpal')
# Sometimes I fix things and theyre only available on the GitHub version
# Right now there are some things that run a lot faster on teh GitHub version
# so I recommend installing that.
# updated 22-12-6 to include loop for PG event adding w/PAMmisc
devtools::install_github('TaikiSan21/PAMpal')
library(PAMpal)

# Start by creating a "PAMpalSettings" object. This keeps track of what data
# you want to process and what processing you want to apply to it.

# Change paths below to your DB and binary folder. Can just be the 
# highest level binary folder for that drift - it will add all files
# within that folder recursively through subfolders.

# This will also ask you to type in some parameters for calculations
# in your console. You can just hit ENTER to accept defaults for all
# of these, they aren't relevant to the GPL calculations only for clicks.

pps <- PAMpalSettings(db = 'path/to/db.sqlite3',
                      binaries = 'path/to/binaryfolder/',
                      # these parameters are only for the click detector - can ignroe
                      sr_hz='auto',
                      filterfrom_khz=0,
                      filterto_khz=NULL,
                      winLen_sec=.0025)

# Now tell it to process your data. Id is optional and serves no function,
# but can be useful to tell data apart at a later point in time. Here
# mode = 'recording' tells it how to organize your data. Most of the time
# we are working with data that have been marked manually into events, 
# so PAMpal wants to organize things into events. mode='db' uses the events
# in the database, and only processes the detectoins you've marked out.
# In this case we just want to process everything, which is what 
# mode='recording' does. It will group them into events by recording file. 

# This might take some time
data <- processPgDetections(pps, mode='recording', id='Humpback007')

# And here's how you can get the detections information out of "data"
# as a dataframe. Time column is "UTC", other columns are stuff it
# measured. 
gplDf <- getGPLData(data)

# Now we can add the wav files to this data. You might get a warning about
# "startSample", its safe to ignore that.
data <- addRecordings(data, folder='path/to/wavfolder')

# that data is stored here as a dataframe. Has "start" & "end" as POSIXct and
# the fulle path to the file as "file"
wavDf <- files(data)$recordings
# add number of detections to this
nDets <- sapply(events(data), nDetections)
nDets <- data.frame(join=names(nDets), nDets=nDets)
wavDf$join <- basename(wavDf$file)
wavDf <- left_join(wavDf, nDets)
wavDf$join <- NULL
wavDf$nDets[is.na(wavDf$nDets)] <- 0

# If we care about assigning some kind of initial label to these
# detections. Otherwise ignore. 
data <- setSpecies(data, method='manual', value='InitialGPL')

# Add events from wavDf loop
for(e in 1:nrow(wavDf)) {
    thisEv <- data[[basename(wavDf$file[e])]]
    # this will get all detector types, if just one type is wanted can
    # be simplified to ex. uids <- unique(getGPLData(thisEv)$UID)
    uids <- unique(unlist(lapply(getDetectorData(thisEv), function(x) {
        if(is.null(x)) return(NULL)
        x$UID
    })))
    
    addPgEvent(db = files(thisEv)$db,
               binary = files(thisEv)$binaries,
               eventType = species(thisEv)$id,
               UIDs = uids,
               type = 'dg',
               start = wavDf$start[e],
               end = wavDf$end[e],
               comment = paste0('Added by PAMpal, event ID: ', id(thisEv)))
}