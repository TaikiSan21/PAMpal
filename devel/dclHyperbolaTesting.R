# sim tdoa / angle

simPath <- function(ship, target, end, n=100, plot=T, changeDepth=FALSE) {
    path <- matrix(c(
        seq(ship[1], end[1], length.out=n),
        seq(ship[2], end[2], length.out=n),
        seq(ship[3], end[3], length.out=n)
    ), ncol=3)
    # print(path[1:3, ])
    arrSep <- .5
    arrPath <- path
    arrPath[, 2] <- arrPath[, 2] - arrSep
    diffs <- numeric(n)
    if(changeDepth) {
        # tarDepths <- c(seq(from=0, to=target[3], length.out=n %/% 2),
        #                seq(from=target[3], length.out = (n %/% 2) + (n %% 2)))
        tarDepths <- c(
            seq(from=0, to=target[3], length.out = n %/% 3),
            rep(target[3], n %/% 3),
            seq(from=target[3], to=0, length.out = (n %/% 3) + (n %% 3))
        )
    } else {
        tarDepths <- rep(target[3], n)
    }
    for(i in 1:nrow(path)) {
        # target[3] <- tarDepths[i]
        dist1 <- sqrt(sum((path[i, ] - c(target[1:2], tarDepths[i]))^2))
        dist2 <- sqrt(sum((arrPath[i, ] - c(target[1:2], tarDepths[i]))^2))
        diffs[i] <- dist2 - dist1
        # diffs[i] <- dist1 - dist2
    }
    a2 <- (diffs/2) ^2
    angle <- sqrt(a2) / (sqrt((arrSep/2)^2 - a2)) * sign(diffs)
    if(plot) {
        plot(x=c(path[, 1], target[1]), y=c(path[, 2], target[2]),
             xlim=c(0, 3e3))
        for(i in 1:n) {
            abline(a=path[i, 2], b=angle[i])
        }
        points(x=target[1], y=target[2], col='red', cex=2)
        points(x=sqrt(target[1]^2+target[3]^2), y=target[2], col='blue', cex=1.5)
    }
    atan(angle) * 180 /pi
}
whale <- c(2e3, 1e3, 0)
ship <- c(0, 0, 0) #0, -.1, 0 array?
end <- c(0, 3e3, 0)

surface <- simPath(ship, whale, end, n=20, changeDepth = T)
deep <- simPath(ship, whale - c(0, 0, -1000), end, n=20, changeDepth = F)
# okay so this definitely only works if depth is fixed, and in this
# case the tarmo loc definitely matches at the point that we correctly
# undo with pythag after figuring out depth. But for very changing depth
# this is not going to happen

# ACTUALLY it might be fine based on sims for a very simple dive profile?
# or its at least not wrong by a lot

# Okay using a pretend dive-stay-surface profile its less clear