context('Test S4 accessors')

test_that('Test AcousticEvent accessors', {
    exPps <- new('PAMpalSettings')
    exPps <- addDatabase(exPps, system.file('extdata', 'Example.sqlite3', package='PAMpal'))
    exPps <- addBinaries(exPps, system.file('extdata', 'Binaries', package='PAMpal'))
    exPps <- addFunction(exPps, roccaWhistleCalcs, module='WhistlesMoans')
    exPps <- addFunction(exPps, standardCepstrumCalcs, module = 'Cepstrum')
    exData <- processPgDetections(exPps, mode='db', id='Example')
    exData <- setSpecies(exData, method='manual', value=letters[1:2])
    exEvent <- exData[[1]]

    expect_equal(id(exEvent), 'Example.OE1')
    dets <- detectors(exEvent)
    # still 3 because we pull UIDs and times for clicks
    expect_equal(length(dets), 3)
    expect_is(dets[[1]], 'data.frame')
    expect_equal(species(exEvent)$id, 'a')
    expect_identical(exEvent[[1]], dets[[1]])
    expect_is(settings(exEvent), 'list')
    expect_is(files(exEvent), 'list')
    expect_equal(files(exEvent)$db, system.file('extdata', 'Example.sqlite3', package='PAMpal'))
    expect_is(ancillary(exEvent), 'list')
})

test_that('Test AcousticStudy accessors', {
    exPps <- new('PAMpalSettings')
    exPps <- addDatabase(exPps, system.file('extdata', 'Example.sqlite3', package='PAMpal'))
    exPps <- addBinaries(exPps, system.file('extdata', 'Binaries', package='PAMpal'))
    exPps <- addFunction(exPps, roccaWhistleCalcs, module='WhistlesMoans')
    exPps <- addFunction(exPps, standardCepstrumCalcs, module = 'Cepstrum')
    exData <- processPgDetections(exPps, mode='db', id='Example')

    expect_identical(exData[[1]], events(exData)[[1]])
    expect_identical(getDetectorData(exData), detectors(exData))
    expect_is(ancillary(exData), 'list')
    expect_is(files(exData), 'list')
    expect_equal(files(exData)$db, system.file('extdata', 'Example.sqlite3', package='PAMpal'))
    expect_identical(pps(exData), exPps)
    expect_is(gps(exData), 'data.frame')
    expect_equal(id(exData), 'Example')
    expect_is(models(exData), 'list')
    expect_is(effort(exData), 'data.frame')
})
