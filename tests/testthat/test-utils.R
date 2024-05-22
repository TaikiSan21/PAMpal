context('Test utils')

test_that('Both ends strtrim', {
    strings <- c('a b', '    a b', '   a   b    ', '    ', 'a  b     ')
    answers <- c('a b', 'a b', 'a   b', '', 'a  b')
    results <- strsplitboth(strings)
    expect_identical(results, answers)
})

test_that('Dets back to study', {
    data(exStudy)
    clicks <- getClickData(exStudy)
    whistles <- getWhistleData(exStudy)
    ceps <- getCepstrumData(exStudy)
    addClicks <- detDataToStudy(exStudy, clicks)
    expect_identical(whistles, getWhistleData(addClicks))
    
    addWhistles <- detDataToStudy(exStudy, whistles)
    expect_equivalent(addWhistles, exStudy)
    addCeps <- detDataToStudy(exStudy, ceps)
    expect_equivalent(addCeps, exStudy)
    clicks$newValue <- 1
    addClicks <- detDataToStudy(exStudy, clicks)
    expect_equivalent(clicks, getClickData(addClicks)[colnames(clicks)])
    
})
