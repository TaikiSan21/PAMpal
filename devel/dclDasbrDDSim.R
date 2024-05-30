# dasbr sim

simDasbr <- function(d1, d2, angle=45, dist=4e3, n=10, plot=FALSE) {

    mid <- -1*mean(c(d1, d2))
    slope <- 1/tan((180-angle)/180*pi)
    xmax <- dist / sqrt(1+slope^2)
    # if(abs(slope) > 1) {
    #     dist <- dist / abs(slope)
    # }
    xs <- seq(from=0, to=xmax, length.out=n)[-1]
    ys <- mid + xs * slope
    if(plot) {
        plot(x=c(0, 0), y=c(-d1, -d2), col='black', xlim=c(0, dist), ylim=c(-dist, 0))
        points(x=xs, y=ys, col='blue')
    }
    t1 <- sqrt(xs^2 + (ys+d2)^2)
    t2 <- sqrt(xs^2 + (ys - d2)^2)
    t2-t1
    list(x=xs, y=ys, td=t2-t1)
}

simDasbr(100, 115, dist=2e3, angle=80, plot=T)

angles <- seq(from=20, to=75, by=5)
simDas <- lapply(angles, function(x) {
    result <- simDasbr(100, 115, x, dist=4e3, n=100)
    result$angle <- x
    result$td <- result$td / 1500
    result
}) %>% bind_rows()

ggplot(simDas) +
    geom_line(aes(x=td*1e3, y=y, col=factor(angle))) +
    ggtitle('Depth vs Time Delay by Angle (90 == Below)') +
    labs(x='Time Delay (ms)', y='Depth (m)')
