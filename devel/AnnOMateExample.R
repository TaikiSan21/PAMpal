#### Annomate Example ####
# Updated 2023-02-06
# requires PAMpal 0.18.0 and PAMmisc 1.11.0

# data must have GPS and species assigned
library(PAMpal)
# data <- readRDS('path/to/AcousticStudy.rds')
specMap <- data.frame(
    old = unique(species(data)),
    new = NA)
specMap$new <- c('write', 'in', 'species', 'reassignments')
anno <- prepAnnotation(data, specMap = specMap, mode='event',
                       # following args are optional, fill in
                       # any that will be the same for all etnries
                       source = 'figshare.com',
                       contact = 'taiki.sakai@noaa.gov')
# or you can add them later
anno$annotator = 'mr taiki'
library(PAMmisc)
# to get figshare data you will need the article id and your personal token
# i recommend storing your token in an external file instead of typing it in here
# also DO NOT commit this file to any repository
figToken <- yaml::read_yaml('secrets.yaml')$fig_token
figId <- 13393360
figshareData <- getFigshareInfo(figToken, figId)

anno <- matchRecordingUrl(anno, figshareData)
# will give you warnings and messages about missing fields
checkAnnotation(anno)
# if you need to fix any manually, write to csv then read back in
annoFile <- 'AnnomateExport.csv'
write.csv(anno, file=annoFile, row.names = FALSE)
# go fix stuff
anno <- read.csv(annoFile, stringsAsFactors = FALSE)
# recheck to see if all seems fine
checkAnnotation(anno)
# add to acoustic study for storage
data <- addAnnotation(data, anno)
# this creates CSV ready for figshare, will also repeat messages from the check
export_annomate(data, file=annoFile)