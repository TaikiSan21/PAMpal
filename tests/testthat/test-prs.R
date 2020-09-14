context('PAMpalSettings object creation and manipulation')

test_that('Create and modify PPS', {
    exPps <- new('PAMpalSettings')
    exPps <- addDatabase(exPps, system.file('extdata', 'Example.sqlite3', package='PAMpal'))
    exPps <- addBinaries(exPps, system.file('extdata', 'Binaries', package='PAMpal'))
    exClick <- function(data) {
        standardClickCalcs(data, filterfrom_khz = 0)
    }
    exPps <- addFunction(exPps, exClick, module = 'ClickDetector')
    exPps <- addFunction(exPps, roccaWhistleCalcs, module='WhistlesMoans')
    exPps <- addFunction(exPps, standardCepstrumCalcs, module = 'Cepstrum')

    expect_is(exPps, 'PAMpalSettings')
    expect_equal(length(exPps@db), 1)
    expect_equal(length(exPps@binaries$list), 3)
    expect_equal(length(exPps@binaries$folder), 1)
    expect_equal(length(exPps@functions), 3)

    exPps <- removeFunction(exPps, 1)
    expect_equal(length(exPps@functions$ClickDetector), 0)
    exPps <- removeDatabase(exPps, 1)
    expect_equal(length(exPps@db), 0)
    exPps <- removeBinaries(exPps, 1)
    expect_equal(length(exPps@binaries$folder), 0)
    expect_equal(length(exPps@binaries$list), 0)

    # calibration
    calFile <- system.file('extdata', 'calibration.csv', package='PAMpal')
    # adding but function doesnt have calibration, hshouldnt fail
    exPps <- addCalibration(exPps, calFile=calFile, all=TRUE, units=3)
    expect_equal(length(exPps@calibration[[1]]), 1)
    # add funciton with cal option, now should work
    calClick <- function(data, calibration=NULL) {
        standardClickCalcs(data, calibration=calibration, filterfrom_khz = 0)
    }
    exPps <- addFunction(exPps, calClick, module = 'ClickDetector')
    exPps <- applyCalibration(exPps, all=TRUE)
    expect_equal(formals(exPps@functions$ClickDetector[[1]])$calibration,
                 calFile)
})
