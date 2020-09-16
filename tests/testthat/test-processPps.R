context('Test processing data with a PPS')

test_that('Test process database', {
    exPps <- new('PAMpalSettings')
    exPps <- addDatabase(exPps, system.file('extdata', 'Example.sqlite3', package='PAMpal'))
    exPps <- addBinaries(exPps, system.file('extdata', 'Binaries', package='PAMpal'))
    exClick <- function(data) {
        standardClickCalcs(data, calibration=NULL, filterfrom_khz = 0)
    }
    exPps <- addFunction(exPps, exClick, module = 'ClickDetector')
    exPps <- addFunction(exPps, roccaWhistleCalcs, module='WhistlesMoans')
    exPps <- addFunction(exPps, standardCepstrumCalcs, module = 'Cepstrum')
    exData <- processPgDetections(exPps, mode='db')

    expect_is(exData, 'AcousticStudy')
    expect_is(exData[1], 'AcousticStudy')
    expect_is(exData[[1]], 'AcousticEvent')
    expect_equal(length(detectors(exData[[1]])), 3)
    # check correct number of dets
    expect_equal(nrow(detectors(exData[[1]])[[1]]), 2)
    expect_equal(nrow(detectors(exData[[1]])[[2]]), 5)
    expect_equal(nrow(detectors(exData[[1]])[[3]]), 7)
    # check no NAs in calcs
    expect_true(!any(
        is.na(detectors(exData[[1]])[[1]]$peak)
    ))
    expect_true(!any(
        is.na(detectors(exData[[1]])[[2]]$ici)
    ))
    expect_true(!any(
        is.na(detectors(exData[[1]])[[3]]$freqBeg)
    ))
})

test_that('Test working with AcousticStudy object', {
    # build basic study object
    exPps <- new('PAMpalSettings')
    exPps <- addDatabase(exPps, system.file('extdata', 'Example.sqlite3', package='PAMpal'))
    exPps <- addBinaries(exPps, system.file('extdata', 'Binaries', package='PAMpal'))
    exClick <- function(data) {
        standardClickCalcs(data, calibration=NULL, filterfrom_khz = 0)
    }
    exPps <- addFunction(exPps, exClick, module = 'ClickDetector')
    exPps <- addFunction(exPps, roccaWhistleCalcs, module='WhistlesMoans')
    exPps <- addFunction(exPps, standardCepstrumCalcs, module = 'Cepstrum')
    exData <- processPgDetections(exPps, mode='db')

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
    exData <- calculateICI(exData)
    expect_true(all(c('Click_Detector_1', 'All') %in% names(ancillary(exData[[1]])$ici)))
    expect_true(!any(
        is.na(ancillary(exData[[1]])$ici[[1]]$ici)
    ))
    expect_true(!any(
        is.na(ancillary(exData[[1]])$ici[[2]]$ici)
    ))
    expect_true(all(c('Click_Detector_1_ici', 'All_ici') %in% names(ancillary(exData[[1]])$measures)))

    # check setSpecies
    exData <- setSpecies(exData, method='pamguard')
    expect_equal(species(exData[[1]])$id, 'Test')
    expect_equal(species(exData[[2]])$id, 'Test')
    # check manual edge cases
    expect_warning(setSpecies(exData, method='manual'), 'Manual mode requires')
    expect_warning(setSpecies(exData, method='manual', value=1:3), 'Length of "value"')
    expect_warning(setSpecies(exData, method='manual', value= data.frame(old=1, new=2),
                              'must contain columns'))
    expect_warning(setSpecies(exData, method='manual',
                              value = data.frame(event = 'a', species=1)),
                   'No match found')
    exData <- setSpecies(exData, method = 'manual', value=letters[1:2])
    expect_equal(species(exData[[1]])$id, 'a')
    expect_equal(species(exData[[2]])$id, 'b')
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
    banterData <- export_banter(exData)
    expect_equal(nrow(banterData$events), 2)
    expect_equal(length(banterData$detectors), 3)
    expect_error(export_banter(exData, dropSpecies = 'b'))
    lessData <- export_banter(exData, dropVars = c('peak'))
    expect_true(!any(
        sapply(lessData$detectors, function(x) 'peak' %in% colnames(x))
    ))

    # test filtering
    filterNone <- filter(exData, VARDNE == 'DNE')
    expect_identical(exData, filterNone)

})

test_that('Test checkStudy test cases', {
    # create example data
    exPps <- new('PAMpalSettings')
    exPps <- addDatabase(exPps, system.file('extdata', 'Example.sqlite3', package='PAMpal'))
    exPps <- addBinaries(exPps, system.file('extdata', 'Binaries', package='PAMpal'))
    exClick <- function(data) {
        standardClickCalcs(data, calibration=NULL, filterfrom_khz = 0)
    }
    exPps <- addFunction(exPps, exClick, module = 'ClickDetector')
    exPps <- addFunction(exPps, roccaWhistleCalcs, module='WhistlesMoans')
    exPps <- addFunction(exPps, standardCepstrumCalcs, module = 'Cepstrum')
    exData <- processPgDetections(exPps, mode='db')
    exData$Example.OE1$Click_Detector_1$peak <- 0
    expect_warning(checkStudy(exData), 'Some clicks had a peak frequency of 0')
})
