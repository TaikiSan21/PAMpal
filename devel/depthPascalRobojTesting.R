es <- readRDS('../Data/RoboJ/Test/eventSummary.rds')
det <- readRDS('../Data/RoboJ/Test/detFunction.rds')

sm3 <- det$SM3M
View(sm3)
dt <- get('dt', envir=environment(sm3$likeFun))

lf <- sm3$likeFun(sm3$Angle)

out <- lf(c(2e3, 1), like=FALSE)

par(mfrow=c(1,1))
# llarr <- makeLLArr(FUN=lf)
# saveRDS(llarr, '../Data/RoboJ/Test/LLArr.rds')
llarr <- readRDS('../Data/RoboJ/Test/LLArr.rds')
q <- quantile(llarr, 0.9)
llarr[llarr < q] <- q
image(x=seq(from=1e3, to=3e3, by=10),
      y=seq(from=-2, to=28, by=.3),
      z=llarr)
points(x=sm3$MaxLikeParam[1], y=sm3$MaxLikeParam[2])



esm3 <- filter(es, es$Recorder == 'SM3M')
colnames(esm3)
esm3 <- rename(esm3, meanAngle = MeanAngle, recorder=Recorder)
detFun <- newEstDetFunction(esm3, subsetCol = NULL, nSim=1e7)

pes <- read.csv('../Data/RoboJ/Test/AllEventClicks_RoboJ_wDrift13K.csv', stringsAsFactors = F)
pes <-readRDS('../Data/RoboJ/Test/eventClicksRC.rds')
pes <- rename(pes, angleCorrected = AngleCorrected, station=Station, species=eventType)
pes <- formatEventSummary(pes, pascal=TRUE)
pes <- filter(pes, station <= 22)
pdf <- newEstDetFunction(pes, doJackknife = F, nSim = 1e7)

comb <- rbind(esm3[c('meanAngle', 'recorder', 'nClicks', 'species')],
              pes[c('meanAngle', 'recorder', 'nClicks', 'species')])
comb$survey <- 'PASCAL'
comb$survey[1:nrow(esm3)] <- 'CCES'
comb %>%
    filter(recorder %in% c('SM3M')) %>%
    ggplot(aes(x=meanAngle, col=survey)) +
    # geom_histogram() +
    geom_density() +
    # facet_wrap(~survey) +
    xlim(0, 90) +
    labs(title='SM3M Angles')

cdf <- newEstDetFunction(comb, doJackknife = FALSE, nSim=1e7)
plotDetFun(c(2836, -1.7), col='black') # running old jay code on CCES, -1.7 is range cap
plotDetFun(detFun$All$maxLikeParam, col='blue', add=TRUE, lwd=2) # my code on CCES, -2 is range cap
plotDetFun(cdf$SM3M$maxLikeParam, col='green', add=TRUE, lwd=2) # combined
plotDetFun(pdf$SM3M$maxLikeParam, col='orange', add=TRUE, lwd=2) # just sm3 from pascal
lf <- detFun$All$likeFun(detFun$All$angle)

llAgain <- makeLLArr(FUN=lf)
saveRDS(llAgain, '../Data/RoboJ/Test/LLAgain.rds')
image(x=seq(from=1e3, to=3e3, by=10), y=seq(from=-2, to=28, by=.3),z=llAgain)
par(mfrow=c(1,2))
q <- quantile(llarr, 0.9)
llarr[llarr < q] <- q
image(llarr)
q <- quantile(llAgain, 0.9)
llAgain[llAgain < q] <- q
image(llAgain)

p1 <- seq(from=1e3, to=3e3, by=50)
p2 <- seq(from=-2, to=5, by=.5)
llDiff <- makeLLArr(p1s = p1,
                    p2s= p2,
                    FUN = lf)
saveRDS(llDiff, file='../Data/RoboJ/Test/llDiff.rds')
llDiff <- readRDS('../Data/RoboJ/Test/llDiff.rds')
q <- quantile(llDiff, 0.8)
llDiff[llDiff < q] <- q
image(x=p1, y=p2, z=-llDiff)

cols <- c('black', 'red', 'blue', 'green')
p1 <- (1:5)*1e3
plotDetFun(c(p1[1], -2))
for(j in seq_along(p1)) {
    plotDetFun(c(p1[j], 1), col=cols[j], add=TRUE)
    for(i in seq(from=-2, to=1, by=.2)) {
        plotDetFun(c(p1[j], i), add=T, col=cols[j])
    }
}

# HP checking
# ec <- readRDS('../Data/RoboJ/Test/eventClicksRC.rds')
# ec$Angle <- ec$Angle + 90
# sspList <- readRDS('../Data/RoboJ/Test/sspList.rds')
data <- readRDS('../Data/RoboJ/Test/RoboJayStudy.rds')
data <- setSpecies(data, method='pamguard')
speciesRename <- data.frame(old = c('BW34-50', 'BW46', 'BW50-75', '?BW', 'Zc', 'Pm', 'SW', '?Pm', 'Oo'),
                            new = c('BW37V', 'MS', 'MS', 'BWunid', 'ZC', 'PM', 'PM', '?PM', 'OO'))
data <- setSpecies(data, method = 'reassign', value=speciesRename)
# Make sure this doesnt have any codes that didnt get converted above
unique(species(data))
data <- filter(data, species != 'Ship')
CCESPat <- '.*Drift-([0-9]{1,2})_.*[\\.OE|\\.DGL][0-9]+$'
ec <- formatClickDf(data, pascal=F, stationPattern = CCESPat)
ec <- filter(ec, species == 'ZC', station != '15')
ec <- splitSSPGroups(ec, splitBy='station', time=3600*24*1)
# mislabeled
ec$station[ec$station == '15'] <- '14'
outPath <- '../Data/RoboJ/Test/'
sspList <- createSSPList(ec, file=file.path(outPath, paste0(as.character(Sys.Date()), '_sspList.rds'))) #may need to add date
sspList <- readRDS('../Data/RoboJ/Test/2022-06-01_sspList.rds')
# ec10 <- doAngleCorrection(ec, sspList, hpDepth=10)
ec100 <- doAngleCorrection(ec, sspList, hpDepth=100)
ec150 <- doAngleCorrection(ec, sspList, hpDepth = 150)
ec200 <- doAngleCorrection(ec, sspList, hpDepth=200)
ecList <- readRDS('ecList.rds')
ec100 <- ecList[[1]]
ec150 <- ecList[[2]]
ec200 <- ecList[[3]]
# ec10$depth <- '10'
ec100$depth <- '100'
ec150$depth <- '150'
ec200$depth <- '200'
allec <- rbind(ec100, ec150, ec200)
ggplot(allec, aes(x=angleCorrected, col=depth)) +
    geom_density()

recorderInfo <- formatDeployDetails('../Data/RoboJ/Test/deployDetails.csv')
recorderInfo <- filter(recorderInfo, Project == 'CCES')
# Split Station 4 - first 20days 2/2, then 2/18
st4new <- recorderInfo[recorderInfo$station == '4', ]
st4new$station <- '4_2'
st4new$dutyCycle <- '2/18'
st4new$deployTime <- st4new$deployTime + 20 * 24 * 3600
recorderInfo <- rbind(recorderInfo, st4new)
recorderInfo$station[recorderInfo$station == '4'] <- '4_1'
recorderInfo$dutyCycle[recorderInfo$station == '4_1'] <- '2/2'
recorderInfo$retrTime[recorderInfo$station == '4_1'] <- recorderInfo$deployTime[recorderInfo$station == '4_2']
# Make changes for eventClicks data. Set all to 4_1, then any after 4_2 deploy time to 4_2
ec$station[ec$station == '4'] <- '4_1'
ec$station[ec$station == '4_1' & ec$UTC >= st4new$deployTime] <- '4_2'

es100 <- formatEventSummary(ec100, recorderInfo = recorderInfo, pascal=FALSE)
es150 <- formatEventSummary(ec150, recorderInfo = recorderInfo, pascal=FALSE)
es200 <- formatEventSummary(ec200, recorderInfo = recorderInfo, pascal=FALSE)
es100$depth <- '100'
es200$depth <- '200'
es150$depth <- '150'
esAll <- bind_rows(es100, es150, es200)
ggplot(esAll, aes(x=meanAngle, col=depth)) +
    geom_density() +
    facet_wrap(~recorder)
df100 <- newEstDetFunction(es100, subsetCol = 'recorder',
                                 doJackknife = TRUE, jk_nSamples = 5, nSim=1e7, progress=TRUE, hpDepth = 100)
df150 <- newEstDetFunction(es150, subsetCol = 'recorder',
                           doJackknife = TRUE, jk_nSamples = 5, nSim=1e7, progress=TRUE, hpDepth = 150)
df200 <- newEstDetFunction(es200, subsetCol = 'recorder',
                           doJackknife = TRUE, jk_nSamples = 5, nSim=1e7, progress=TRUE, hpDepth = 200)

sapply(df100, function(x) x$maxLikeParam)
sapply(df100, function(x) x$maxL_EDR)
sapply(df150, function(x) x$maxLikeParam)
sapply(df150, function(x) x$maxL_EDR)
sapply(df200, function(x) x$maxLikeParam)
sapply(df200, function(x) x$maxL_EDR)
dfList <- readRDS('dfList.rds')
df100 <- dfList[[1]]
df150 <- dfList[[2]]
df200 <- dfList[[3]]
plotDetFun(df100$SM3M$maxLikeParam, title='SM3M 100/150/200m')
plotDetFun(df150$SM3M$maxLikeParam, col='blue', add=T)
plotDetFun(df200$SM3M$maxLikeParam, col='red', add=T)

plotDetFun(df100$ST4300HF$maxLikeParam, title='ST4300 100/150/200m')
plotDetFun(df150$ST4300HF$maxLikeParam, col='blue', add=T)
plotDetFun(df200$ST4300HF$maxLikeParam, col='red', add=T)

angle100 <- get('dt', envir = environment(df100$SM3M$likeFun))
angle150 <- get('dt', envir = environment(df150$SM3M$likeFun))
angle200 <- get('dt', envir = environment(df200$SM3M$likeFun))
angle100$iAngle <- as.numeric(angle100$iAngle)
angle150$iAngle <- as.numeric(angle150$iAngle)
angle200$iAngle <- as.numeric(angle200$iAngle)
angle100$depth <- '100'
angle150$depth <- '150'
angle200$depth <- '200'
ang <- rbind(angle100, angle150, angle200)
ggplot(ang, aes(x=iAngle, col=depth)) +
    geom_density(adjust=1.5)
d100 <- get('depth', envir = environment(df100$SM3M$likeFun))
d150 <- get('depth', envir = environment(df150$SM3M$likeFun))
summary(d100)
summary(d150)
