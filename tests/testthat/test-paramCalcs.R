context('Test parameters calculations.')

test_that('Test click calcs', {
    data(testClick)
    clickData <- standardClickCalcs(testClick)
    expect_true(all(c('Channel', 'peak', 'peak2', 'Q_10dB') %in% colnames(clickData)))
    expect_is(clickData, 'data.frame')
    expect_equal(ncol(clickData), 23)
    expect_equal(nrow(clickData), 2)
    expect_true(!any(is.na(clickData)))
    expect_true(clickChecker(standardClickCalcs))
})

test_that('Test whistle calcs', {
    data(testWhistle)
    whistleData <- roccaWhistleCalcs(testWhistle)
    expect_true(whistleChecker(roccaWhistleCalcs))
    expect_true(all(c('freqBeg', 'freqPosSlopeMean', 'freqCofm') %in%
                    names(whistleData)))
    expect_true(!any(is.na(whistleData)))

})

test_that('Test cepstrum calcs', {
    data(testCeps)
    cepsData <- standardCepstrumCalcs(testCeps)
    expect_true(cepstrumChecker(standardCepstrumCalcs))
    expect_true(all(c('ici', 'duration', 'iciSlope') %in% names(cepsData)))
    expect_true(!any(is.na(cepsData)))
})
