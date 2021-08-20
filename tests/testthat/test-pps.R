context('PAMpalSettings object creation and manipulation')

test_that('Create and modify PPS', {
    exPps <- new('PAMpalSettings')
    db <- system.file('extdata', 'Example.sqlite3', package='PAMpal')
    exPps <- addDatabase(exPps, db, verbose=FALSE)
    bin <- system.file('extdata', 'Binaries', package='PAMpal')
    exPps <- addBinaries(exPps, bin, verbose=FALSE)
    exClick <- function(data) {
        standardClickCalcs(data, filterfrom_khz = 0)
    }
    exPps <- addFunction(exPps, exClick, module = 'ClickDetector', verbose=FALSE)
    exPps <- addFunction(exPps, roccaWhistleCalcs, module='WhistlesMoans', verbose=FALSE)
    exPps <- addFunction(exPps, standardCepstrumCalcs, module = 'Cepstrum', verbose=FALSE)

    expect_is(exPps, 'PAMpalSettings')
    expect_equal(length(exPps@db), 1)
    expect_equal(length(exPps@binaries$list), 3)
    expect_equal(length(exPps@binaries$folder), 1)
    expect_equal(length(exPps@functions), 3)

    # test removal
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
    exPps <- addFunction(exPps, calClick, module = 'ClickDetector', verbose=FALSE)
    exPps <- applyCalibration(exPps, all=TRUE)
    expect_equal(formals(exPps@functions$ClickDetector[[1]])$calibration,
                 calFile)
    
    # settings
    xmlFile <- system.file('extdata', 'Example.xml', package='PAMpal')
    xmlList <- loadPamguardXML(xmlFile)
    xmlAdd <- addSettings(exPps, xmlFile, type='xml', verbose=FALSE)
    listAdd <- addSettings(exPps, xmlList, type='list', verbose=FALSE)
    xmlPps <- PAMpalSettings(db, bin, sr_hz='auto', filterfrom_khz=5, filterto_khz=NULL, winLen_sec=.0025, settings=xmlFile, verbose=FALSE)
    expect_identical(xmlPps@settings[1:2], xmlAdd@settings[1:2], listAdd@settings[1:2])
    expect_identical(xmlList[1:2], xmlAdd@settings[1:2])
    expect_true(all(c('sources', 'detectors', 'raw') %in% names(xmlAdd@settings)))
    expect_equal(xmlList$detectors[['Click_Detector']]$sr, 192000)
})
