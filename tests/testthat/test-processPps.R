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
    exData <- processPgDetections(exPps, mode='db', id='Example', progress=FALSE)

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
    expect_warning(processPgDetections(exPps, mode='db', id='Test', progress=FALSE),
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
    exTime <- processPgDetections(exPps, mode='time', grouping=grp, id='Time', progress=FALSE)
    grpChar <- data.frame(start = c('2018-03-20 15:25'),
                      end = c('2018-03-20 15:25:11'),
                      id = 'TimeExample')
    exChar <- processPgDetections(exPps, mode='time', grouping=grpChar, id='Time', format='%Y-%m-%d %H:%M:%S', progress=FALSE)
    expect_identical(exChar, exTime)
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
