context('Test S4 accessors')

test_that('Test AcousticEvent accessors', {
    data(exStudy)
    exStudy <- setSpecies(exStudy, method='manual', value=letters[1:2])
    exEvent <- exStudy[[1]]

    expect_equal(id(exEvent), 'Example.OE1')
    dets <- detectors(exEvent)
    # still 3 because we pull UIDs and times for clicks
    expect_equal(length(dets), 3)
    detectors(exEvent) <- dets[1:2]
    expect_equal(length(detectors(exEvent)), 2)
    expect_is(dets[[1]], 'data.frame')
    expect_equal(species(exEvent)$id, 'a')
    species(exEvent)$id <- 'new'
    expect_equal(species(exEvent)$id, 'new')
    expect_identical(exEvent[[1]], dets[[1]])
    expect_is(settings(exEvent), 'list')
    expect_is(files(exEvent), 'list')
    expect_is(ancillary(exEvent), 'list')
    ancillary(exEvent) <- list(test='new')
    expect_identical(ancillary(exEvent), list(test='new'))
})

test_that('Test AcousticStudy accessors', {
    data(exStudy)
    expect_identical(exStudy[[1]], events(exStudy)[[1]])
    expect_identical(getDetectorData(exStudy), detectors(exStudy))
    expect_is(ancillary(exStudy), 'list')
    ancillary(exStudy) <- list(test='new')
    expect_identical(ancillary(exStudy), list(test='new'))
    expect_is(files(exStudy), 'list')
    expect_is(gps(exStudy), 'data.frame')
    expect_equal(id(exStudy), 'ExampleData_10-12-2020')
    expect_is(models(exStudy), 'list')
    expect_is(effort(exStudy), 'data.frame')
})
