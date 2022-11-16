det <-  readRDS(file.choose())

lf <- det$All$likeFun

hm <- lf(environment(det$All[["likeFun"]])[["Angle"]])

outs <- hm(c(2e3,1), F)
list(Angle=1:91, Density=tf$AngleDensity, nSim=nSim[i])
library(ggplot2)

plotters <- data.frame(Angle=1:91, AngleDensity=outs$AngleDensity)
ggplot(plotters, aes(x=Angle, y=AngleDensity)) + geom_line()

angles <- data.frame(Angle = factor(ceiling(outs$Angle)), SRange=outs$SRange)
ggplot(angles, aes(x=SRange, y=Angle, fill=..x..)) +
    geom_density_ridges_gradient() +
    xlim(0, 4e3) +
    scale_fill_viridis()


dfOut <- denOut
for(i in seq_along(nSim)) {
    simFun <- createSimLikeFun(DepthDistr='log-normal', Model='C_HN', nSim=nSim[i])(angle)
    tf <- simFun(c(1e3, 1), F)
    dfOut[[i]] <- data.frame(Angle=factor(ceiling(tf$Angle)), SRange=tf$SRange, nSim=nSim[i])
    denOut[[i]] <- list(Angle=1:91, Density=tf$AngleDensity, nSim=nSim[i])
}
denOut <- bind_rows(denOut)
dfOut <- bind_rows(dfOut)
gDen <- ggplot(denOut, aes(x=Angle, y=Density, col=factor(nSim))) + geom_line(size=1)
# gRidge <- ggplot(dfOut) +
#     geom_density_ridges_gradient(aes(x=SRange, y=Angle, fill=..x..)) +
#     scale_fill_viridis() +
#     xlim(0, 4e3) +
#     facet_wrap(~nSi

wat <- newEstDetFunction(es, doJackknife = FALSE, nSim=1e7, subsetCol = 'Recorder')
plotDetFun(wat$All$maxLikeParam)
plotDetFun(wat$ST4300HF$maxLikeParam, col='blue', add=TRUE)
plotDetFun(wat$SM3M$maxLikeParam, col='red', add=TRUE)

angles <- bind_rows(list(Recorder='SM3M', Angle=wat$SM3M$angle),
                list(Recorder='ST4300', Angle=wat$ST4300HF$angle))
ggplot(angles, aes(x=Angle, col=Recorder)) + geom_density() +
    geom_vline(xintercept = 75)
