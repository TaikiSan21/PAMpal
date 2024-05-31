# dive depth error sim testing

calculateDepth(arrDepth=5.5, slantRange=2e3, delayTime = .002, soundSpeed = 1500)

estDepthError <- function(arrDepth=5.5, slantRange=2e3, delayTime=.002, soundSpeed=1500,
                          arrErr=1, slantErr=.1, delayErr=.1, ssErr = .05,
                          simFun=c('unif', 'norm'), n=1e3) {
    switch(match.arg(simFun),
           'unif' = {
               arrSims <- arrDepth + runif(n, min=-arrErr, max=arrErr)
               slantSims <- slantRange * (1 + runif(n, min=-slantErr, max=slantErr))
               delaySims <- delayTime * (1 + runif(n, min=-delayErr, max=delayErr))
               ssSims <- soundSpeed * (1 + runif(n, min=-ssErr, max=ssErr))
           },
           'norm' = {
               arrSims <- rnorm(n, mean=arrDepth, sd=arrErr)
               slantSims <- rnorm(n, mean=1, sd=slantErr) * slantRange
               delaySims <- rnorm(n, mean=1, sd=delayErr) * delayTime
               ssSims <- rnorm(n, mean=1, sd=ssErr) * soundSpeed
           }
    )
    
    
    depths <- numeric(length=n)
    for(i in seq_along(depths)) {
        depths[i] <- calculateDepth(arrDepth=arrSims[i],
                                    slantRange = slantSims[i],
                                    delayTime = delaySims[i],
                                    soundSpeed = ssSims[i])
    }
    c(calculateDepth(arrDepth, slantRange, delayTime, soundSpeed), depths)
}


testDepths <- estDepthError(arrDepth=8, arrErr=1, delayTime=.004, n=1000, simFun='norm')
quantile(testDepths, c(.025, .975))
hist(testDepths, breaks=20)
lines(x=rep(testDepths[1], 2), y=c(0, 1e4))

hm <- getClickData(data)

checkErrSource <- function(arrDepth, slantRange, delayTime, soundSpeed,
                           arrErr, slantErr, delayErr, ssErr,
                           sumFun=c('unif', 'norm'), n=1e3) {
    arrOnly <- estDepthError(arrDepth, slantRange, delayTime, soundSpeed,
                             arrErr=arrErr, slantErr=0, delayErr=0, ssErr=0)
    slantOnly <- estDepthError(arrDepth, slantRange, delayTime, soundSpeed,
                             arrErr=0, slantErr=slantErr, delayErr=0, ssErr=0)
    delayOnly <- estDepthError(arrDepth, slantRange, delayTime, soundSpeed,
                             arrErr=0, slantErr=0, delayErr=delayErr, ssErr=0)
    ssOnly <- estDepthError(arrDepth, slantRange, delayTime, soundSpeed,
                             arrErr=0, slantErr=0, delayErr=0, ssErr=ssErr)
    print(quantile((arrOnly-mean(arrOnly))/mean(arrOnly), c(.025, .975)))
    print(quantile(slantOnly, c(.025, .975)))
    print(quantile(delayOnly, c(.025, .975)))
    print(quantile(ssOnly, c(.025, .975)))
}
invisible(checkErrSource(arrDepth=8, slantRange=2e3, delayTime=.004, soundSpeed=1500,
                         arrErr=1, slantErr=.1, delayErr=.1, ssErr=.1))
