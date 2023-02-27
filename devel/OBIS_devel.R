source('./devel/OBIS_Functions.R')
data <- readRDS('../../beaker-wigner-class/rds/ccesStudy.rds')
table(species(data))
# must have "ours" for every species in study we want to export record of
# lookup aphia id here https://www.marinespecies.org/aphia.php?p=search
# coord might be different for each species, roughly max detection distance
speciesMap <- tribble(
    ~ours, ~aphia, ~coordinateUncertaintyInMeters,
    'BB', 242608, 4000,
    'BW', 136986, 4000,
    'BW37V', 136986, 4000,
    'BW43', 136986, 4000,
    'BWC', 136986, 4000,
    'MS', 254991, 4000,
    'ZC',137127, 4000
)
# some fields cannot be filled with PAMpal data, so list their values here to have
# those values filled in for every event/occurrence they are needed for
otherFields <- list(
    basisOfRecord = 'MachineObservation',
    institutionCode = 'SWFSC',
    georeferenceProtocol = 'SPOT GPS',
    samplingProtocol = 'https://doi.org/10.6084/m9.figshare.19358036',
    countryCode = 'US',
    geodeticDatum = 'WGS84')
# dbPattern tells how to convert basename(db) into the event ID using
# gsub(dbPattern[1], dbPattern[2], basename(db))
obis <- export_obis(data, 
                    speciesMap, 
                    dbPattern=c('.*(Drift-[0-9]{1,2}).*', 'CCES_\\1'), 
                    otherFields)
# "export_obis" outputs a list of two dataframes, "event" for event data
# and "occurrence" for occurrence data

library(obistools)
# this gives an error message about parentEventID, I dont fully udnerstand
# where this column is supposed to be
obistools::check_eventids(obis$event)
obistools::check_eventdate(obis$event)
obistools::check_eventdate(obis$occurrence)
obistools::check_fields(obis$occurrence)
write.csv(obis$occur, file='../Data/OBIS_Workshop/CCESTest_Occur.csv', row.names = FALSE)
write.csv(obis$event, file='../Data/OBIS_Workshop/CCESTest_Event.csv', row.names = FALSE)
