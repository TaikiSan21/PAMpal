context('Test processing data with a PPS')

test_that('Test process database', {
    exPps <- new('PAMpalSettings')
    exPps <- addDatabase(exPps, system.file('extdata', 'Example.sqlite3', package='PAMpal'), verbose=FALSE)
    exPps <- addBinaries(exPps, system.file('extdata', 'Binaries', package='PAMpal'), verbose=FALSE)
    exClick <- function(data) {
        standardClickCalcs(data, calibration=NULL, filterfrom_khz = 0)
    }
    exPps <- addFunction(exPps, exClick, module = 'ClickDetector', verbose=FALSE)
    exPps <- addFunction(exPps, roccaWhistleCalcs, module='WhistlesMoans', verbose=FALSE)
    exPps <- addFunction(exPps, standardCepstrumCalcs, module = 'Cepstrum', verbose=FALSE)
    
    setArgPps <- removeFunction(exPps, 1)
    setArgPps <- addFunction(setArgPps, standardClickCalcs, sr_hz='auto', filterfrom_khz=0, filterto_khz=NULL,
                             winLen_sec=.0025, module='ClickDetector', verbose=FALSE)
    
    exData <- processPgDetections(exPps, mode='db', id='Example', progress=FALSE, verbose=FALSE)
    setArgData <- processPgDetections(setArgPps, mode='db', id='Example', progress=FALSE, verbose=FALSE)
    
    expect_identical(events(exData), events(setArgData))
    
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
    exPps <- removeFunction(exPps, 1)
    expect_warning(processPgDetections(exPps, mode='db', id='Test', progress=FALSE, verbose=FALSE),
                   'No functions for processing Module Type: ClickDetector')
})

test_that('Test process time', {
    exPps <- new('PAMpalSettings')
    exPps <- addDatabase(exPps, system.file('extdata', 'Example.sqlite3', package='PAMpal'), verbose=FALSE)
    exPps <- addBinaries(exPps, system.file('extdata', 'Binaries', package='PAMpal'), verbose=FALSE)
    exClick <- function(data) {
        standardClickCalcs(data, calibration=NULL, filterfrom_khz = 0)
    }
    exPps <- addFunction(exPps, exClick, module = 'ClickDetector', verbose=FALSE)
    exPps <- addFunction(exPps, roccaWhistleCalcs, module='WhistlesMoans', verbose=FALSE)
    exPps <- addFunction(exPps, standardCepstrumCalcs, module = 'Cepstrum', verbose=FALSE)
    grp <- data.frame(start = as.POSIXct('2018-03-20 15:25:00', tz='UTC'),
                      end = as.POSIXct('2018-03-20 15:25:11', tz='UTC'),
                      id = 'TimeExample')
    exTime <- processPgDetections(exPps, mode='time', grouping=grp, id='Time', progress=FALSE, verbose=FALSE)
    grpChar <- data.frame(start = c('2018-03-20 15:25'),
                      end = c('2018-03-20 15:25:11'),
                      id = 'TimeExample')
    exChar <- processPgDetections(exPps, mode='time', grouping=grpChar, id='Time', format='%Y-%m-%d %H:%M:%S', 
                                  progress=FALSE, verbose=FALSE)
    expect_identical(events(exChar), events(exTime))
    dets <- getDetectorData(exTime)
    times <- do.call(rbind, lapply(dets, function(x) {
        x[, c('UTC', 'UID')]
    }))
    expect_true(all(times$UTC <= grp$end))
    expect_true(all(times$UTC >= grp$start))
    expect_equal(length(events(exTime)), 1)
    expect_equal(id(exTime[[1]]), grp$id)
    expect_equal(id(exTime), 'Time')
})

test_that('Test process recording', {
    exPps <- new('PAMpalSettings')
    exPps <- addDatabase(exPps, system.file('extdata', 'Example.sqlite3', package='PAMpal'), verbose=FALSE)
    exPps <- addBinaries(exPps, system.file('extdata', 'Binaries', package='PAMpal'), verbose=FALSE)
    exClick <- function(data) {
        standardClickCalcs(data, calibration=NULL, filterfrom_khz = 0)
    }
    exPps <- addFunction(exPps, exClick, module = 'ClickDetector', verbose=FALSE)
    exPps <- addFunction(exPps, roccaWhistleCalcs, module='WhistlesMoans', verbose=FALSE)
    exPps <- addFunction(exPps, standardCepstrumCalcs, module = 'Cepstrum', verbose=FALSE)
    
    exRec <- processPgDetections(exPps, mode='recording', id='rec', progress=FALSE, verbose = FALSE)
    allBin <- do.call(rbind, lapply(exPps@binaries$list, function(x) {
        data.frame(loadPamguardBinaryFile(x, convertDate=TRUE, skipLarge=TRUE))[c('UID', 'date')]
    }))
    allDet <- do.call(rbind, lapply(getDetectorData(exRec), function(x) {
        x[c('UTC', 'UID')]
    }))
    expect_true(all(allBin$UID %in% allDet$UID))
})

test_that('Test process fixed time', {
    exPps <- new('PAMpalSettings')
    exPps <- addDatabase(exPps, system.file('extdata', 'Example.sqlite3', package='PAMpal'), verbose=FALSE)
    exPps <- addBinaries(exPps, system.file('extdata', 'Binaries', package='PAMpal'), verbose=FALSE)
    exClick <- function(data) {
        standardClickCalcs(data, calibration=NULL, filterfrom_khz = 0)
    }
    exPps <- addFunction(exPps, exClick, module = 'ClickDetector', verbose=FALSE)
    exPps <- addFunction(exPps, roccaWhistleCalcs, module='WhistlesMoans', verbose=FALSE)
    exPps <- addFunction(exPps, standardCepstrumCalcs, module = 'Cepstrum', verbose=FALSE)
    exNumeric <- processPgDetections(exPps, mode='fixed', id='fixNumeric', grouping=60, progress=FALSE, verbose=FALSE)
    exMin <- processPgDetections(exPps, mode='fixed', id='fixMin', grouping='1min', progress=FALSE, verbose=FALSE)
    expect_identical(exNumeric[[1]], exMin[[1]])
    expect_identical(ancillary(exNumeric)$grouping,
                     ancillary(exMin)$grouping)
})