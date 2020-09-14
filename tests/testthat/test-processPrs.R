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

