context('Test working with AcousticStudy object')

test_that('Test working with AcousticStudy object', {
    # build basic study object
    exPps <- new('PAMpalSettings')
    exPps <- addDatabase(exPps, system.file('extdata', 'Example.sqlite3', package='PAMpal'), verbose=FALSE)
    exPps <- addBinaries(exPps, system.file('extdata', 'Binaries', package='PAMpal'), verbose=FALSE)
    exClick <- function(data) {
        standardClickCalcs(data, calibration=NULL, filterfrom_khz = 0)
    }
    exPps <- addFunction(exPps, exClick, module = 'ClickDetector', verbose=FALSE)
    exPps <- addFunction(exPps, roccaWhistleCalcs, module='WhistlesMoans', verbose=FALSE)
    exPps <- addFunction(exPps, standardCepstrumCalcs, module = 'Cepstrum', verbose=FALSE)
    exData <- processPgDetections(exPps, mode='db', id='Example', progress = FALSE, verbose = FALSE)

    # check adding gps
    exData <- addGps(exData)
    expect_equal(nrow(gps(exData)), 200)
    expect_true(!any(
        is.na(gps(exData)[['Latitude']])
    ))
    expect_true(all(c('UTC', 'Latitude', 'Longitude') %in% colnames(exData[[1]][[1]])))
    expect_true(!any(
        is.na(exData[[1]][[1]][['Latitude']])
    ))
    expect_true(!any(
        is.na(exData[[1]][[1]][['Longitude']])
    ))
    expect_true(all(c('UTC', 'Latitude', 'Longitude') %in% colnames(exData[[1]][[2]])))
    expect_true(!any(
        is.na(exData[[1]][[2]][['Latitude']])
    ))
    expect_true(!any(
        is.na(exData[[1]][[2]][['Longitude']])
    ))
    expect_true(all(c('UTC', 'Latitude', 'Longitude') %in% colnames(exData[[1]][[3]])))
    expect_true(!any(
        is.na(exData[[1]][[3]][['Latitude']])
    ))
    expect_true(!any(
        is.na(exData[[1]][[3]][['Longitude']])
    ))

    # check ici
    expect_warning(iciData <- getICI(exData), 'No ICI data')
    expect_null(iciData)
    exData <- calculateICI(exData)
    expect_true(all(c('Click_Detector_1', 'All') %in% names(ancillary(exData[[1]])$ici)))
    expect_true(!any(
        is.na(ancillary(exData[[1]])$ici[[1]]$ici)
    ))
    expect_true(!any(
        is.na(ancillary(exData[[1]])$ici[[2]]$ici)
    ))
    expect_true(all(c('Click_Detector_1_ici', 'All_ici') %in% names(ancillary(exData[[1]])$measures)))
    expect_true(all(c('Click_Detector_1_ici', 'All_ici') %in% names(getClickData(exData[[1]]))))
    iciData <- getICI(exData, 'data')
    expect_true(all(c('Click_Detector_1', 'All') %in% names(iciData[[1]])))
    expect_identical(names(iciData), names(events(exData)))
    iciData <- getICI(exData, 'value')
    expect_true(all(c('Click_Detector_1_ici', 'All_ici') %in% names(iciData[[1]])))

    # check setSpecies
    exData <- setSpecies(exData, method='pamguard')
    expect_equal(species(exData[[1]])$id, 'Test')
    expect_equal(species(exData[[2]])$id, 'Test')
    # check manual edge cases
    expect_warning(setSpecies(exData, method='manual'), 'Manual mode requires')
    expect_warning(setSpecies(exData, method='manual', value=1:4), 'Length of "value"')
    expect_warning(setSpecies(exData, method='manual', value= data.frame(old=1, new=2),
                              'must contain columns'))
    expect_message(setSpecies(exData, method='manual',
                              value = data.frame(event = 'a', species=1)),
                   'No match found')
    exData <- setSpecies(exData, method = 'manual', value=letters[1:3])
    expect_equal(species(exData[[1]])$id, 'a')
    expect_equal(species(exData[[2]])$id, 'b')
    expect_equal(species(exData[[3]])$id, 'c')
    exData <- setSpecies(exData, method='manual',
                         value = data.frame(event='Example.OE1', species = 'c'))
    expect_equal(species(exData[[1]])$id, 'c')
    # check reassign edge cases
    expect_warning(setSpecies(exData, method='reassign'), 'mode requires a "value"')
    expect_warning(setSpecies(exData, method='reassign', value=data.frame(x=1, y=2)),
                   'must have columns')
    exData <- setSpecies(exData, method='reassign',
                         value= data.frame(old='c', new='b'))
    expect_equal(species(exData[[1]])$id, 'b')
    # test banter export
    banterData <- export_banter(exData, verbose=FALSE)
    expect_equal(nrow(banterData$events), 3)
    expect_equal(length(banterData$detectors), 3)
    expect_warning(export_banter(exData, dropSpecies = c('b', 'c'), verbose=FALSE))
    lessData <- export_banter(exData, dropVars = c('peak'), verbose=FALSE)
    expect_true(!any(
        sapply(lessData$detectors, function(x) 'peak' %in% colnames(x))
    ))

    # test add recordings
    recs <- system.file('extdata', 'Recordings', package='PAMpal')
    exData <- addRecordings(exData, folder = recs, log=FALSE, progress=FALSE)
    expect_identical(normalizePath(files(exData)$recordings$file),
                     normalizePath(list.files(recs, full.names = TRUE)))
    expect_warning(warnRec <- addRecordings(exData, folder = 'DNE', log=FALSE, progress=FALSE))

    # test warning access from recorder warning
    warns <- getWarnings(warnRec)
    expect_is(warns, 'data.frame')
    expect_true('Provided folder DNE does not exist.' %in% warns$message)

})

test_that('Test filter', {
    data(exStudy)
    # test filtering
    filterNone <- filter(exStudy, VARDNE == 'DNE')
    expect_identical(events(exStudy), events(filterNone))
    exStudy <- setSpecies(exStudy, method='manual', value=letters[1:2])
    spFilter <- filter(exStudy, species == 'a')
    expect_equal(length(events(spFilter)), 1)
    expect_equal(species(spFilter[[1]])$id, 'a')
    spFilter <- filter(exStudy, species %in% letters[1:3])
    expect_identical(events(spFilter), events(exStudy))
    peakFilter <- filter(exStudy, peak < 20)
    expect_true(all(detectors(peakFilter)$click$peak < 20))
    peakFilter <- filter(exStudy, peak < 2000)
    events(peakFilter) <- lapply(events(peakFilter), function(x) {
        detectors(x) <- lapply(detectors(x), function(y) {
            row.names(y) <- NULL
            y
        })
        x
    })
    events(exStudy) <- lapply(events(exStudy), function(x) {
        detectors(x) <- lapply(detectors(x), function(y) {
            row.names(y) <- NULL
            y
        })
        x
    })
    expect_identical(events(peakFilter), events(exStudy))

    dbFilter <- filter(exStudy, database == files(exStudy)$db)
    expect_identical(events(exStudy), events(dbFilter))
    dbNone <- filter(exStudy, database == 'NODB.sqlite3')
    expect_equal(length(events(dbNone)), 0)
})
test_that('Test checkStudy test cases', {
    # create example data
    data(exStudy)
    expect_warning(checkStudy(exStudy, maxLength = 1),
                   'Found 2 events longer than 1 seconds')
    expect_warning(checkStudy(exStudy, maxSep = .1),
                   'Found 2 events with detections more than 0.1')
    exStudy$Example.OE1$Click_Detector_1$peak <- 0
    expect_warning(checkStudy(exStudy), 'Some clicks had a peak frequency of 0')
})

test_that('Test getBinaryData', {
    data(exStudy)
    binFolder <- system.file('extdata', 'Binaries', package='PAMpal')
    exStudy <- updateFiles(exStudy, bin=binFolder, db=NA, verbose=FALSE)
    bin <- getBinaryData(exStudy, UID = 8000003)
    expect_equal(names(bin), '8000003')
    expect_true(all(c('wave', 'sr', 'minFreq') %in% names(bin[[1]])))
    expect_null(expect_warning(getBinaryData(exStudy, UID = 1)))
})

test_that('Test getDetectorData', {
    data(exStudy)
    dets <- getDetectorData(exStudy)
    expect_true(all(c('click', 'whistle', 'cepstrum') %in% names(dets)))
    expect_is(dets, 'list')
    expect_is(dets[[1]], 'data.frame')
    # works same on events and studies
    expect_identical(getDetectorData(exStudy[1]),
                     getDetectorData(exStudy[[1]]))
    expect_identical(dets$click, getClickData(exStudy))
    expect_identical(dets$whistle, getWhistleData(exStudy))
    expect_identical(dets$cepstrum, getCepstrumData(exStudy))
})

test_that('Test updateFiles', {
    data(exStudy)
    # corrupting filepaths
    files(exStudy)$db <- substr(files(exStudy)$db, start=5, stop=10e3)
    files(exStudy)$binaries <- substr(files(exStudy)$binaries, start=5, stop=10e3)
    files(exStudy[[1]])$db <- substr(files(exStudy[[1]])$db, start=5, stop=10e3)
    files(exStudy[[1]])$binaries <- substr(files(exStudy[[1]])$binaries, start=5, stop=10e3)
    db <- system.file('extdata', 'Example.sqlite3', package='PAMpal')
    bin <- system.file('extdata', 'Binaries', package='PAMpal')
    expect_true(!any(file.exists(files(exStudy)$db,
                                 files(exStudy)$binaries,
                                 # files(exStudy)$recordings$file,
                                 files(exStudy[[1]])$db,
                                 files(exStudy[[1]])$binaries)))
    exStudy <- updateFiles(exStudy, db=db, bin=bin, verbose=FALSE)
    # exStudy <- updateFiles(exStudy, db=db, bin=bin, recording = recs, verbose=FALSE)
    expect_true(all(file.exists(files(exStudy)$db,
                                files(exStudy)$binaries,
                                # files(exStudy)$recordings$file,
                                files(exStudy[[1]])$db,
                                files(exStudy[[1]])$binaries)))
    recs <- system.file('extdata', 'Recordings', package='PAMpal')
    exStudy <- addRecordings(exStudy, folder =recs, log=FALSE, progress=FALSE)
    files(exStudy)$recordings$file <- substr(files(exStudy)$recordings$file, start=5, stop=10e3)
    expect_true(!any(file.exists(files(exStudy)$recordings$file)))
    exStudy <- updateFiles(exStudy, recording=recs, verbose=FALSE)
    expect_true(all(file.exists(files(exStudy)$recordings$file)))
})
