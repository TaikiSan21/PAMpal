source('addAnnotation.R')
# this will ask questions about whether all same/some could be different/not used at all
# if you select all same, you will be prompted to type in. So have the URL for the dataset
# ready for the recording_info_url and such.
# FOR THE RECORDING URL select the not used option, we'll deal with these manually later
# it will also ask you to provide full names for any species IDs it finds
anno <- prepAnnotation(myData)
# This should run without any help unless you selected the "might be different" option above,
# in which case you will be prompted to enter that info for each event
annoTest <- addAnnotation(myData, anno)

# getting figshare URLs to deal with those automatically
install.packages('rfigshare')
library(rfigshare)
# This step will ask you to login through the website
details <- fs_details(13393360)
figdf <- bind_rows(lapply(details$files, function(x) {
    x[c('download_url', 'name')]
}))
figdf$event <- gsub('Event_(.*)CH[1-6].*', '\\1', figdf$name)

# attaching the URLs from figshare to the proper event
for(e in names(events(annoTest))) {
    if(e %in% figdf$event) {
        ancillary(annoTest[[e]])$annotation$recording_url <-
            figdf$download_url[figdf$event == e]
    }
}

# this creates a dataframe that annomate needs to upload. if you want to
# create the csv now, then do write.csv(annoCsv, file='annomate.csv')
# after this, but since we are still waiting on ed i decided not to create the file
# yet

# you will get warnings for any events that are missing any mandatory fields, and these
# get dropped. Currently one event gets dropped because the recording file for it is not
# on the figshare
annoCsv <- export_annomate(annoTest)
