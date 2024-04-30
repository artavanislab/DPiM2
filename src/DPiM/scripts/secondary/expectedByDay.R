{
    byDay <- read.table("~/DPiM/data/inclusionByDate.dpim3.7-28.tab",
                        header = T)
    p <- sum(byDay$bad) / sum(byDay$bad + byDay$good) # likelihood that any
    ## given run is bad
    mx <- max(byDay$bad)
    tot <- byDay$bad + byDay$good;
    x <- factor(byDay$bad, levels=c(0:mx))
    obser <- table(x)
    #chiTest <- chisq.test(expec, obser)

    N <- length(tot)
    expec <- rep(0, mx+1)
    for (i in 1:N) {
        expec <- expec + dbinom(c(0:mx), tot[i], p);
    }

    fscal <- 1.4
    png(file="~/public_html/DPiM/badDaysHistogram.7-28.png", height=500,
        width=600);
    par(oma=c(0,1,0,0))
    plot(c(0:mx), expec, ylim=c(0, max(expec[1], obser[[1]] ) ),
         xlab="Number of bad runs", ylab="Number of days", 
         cex.lab=1.1*fscal, cex.axis= fscal, lwd=2, col="gray40")
    points(table(byDay$bad), lwd=2)
    legend("topright", legend=c("Observed", "Uniform expectation"),
           pch=c(26, 1),
           lty=c(1, 0), bty="n", lwd=2, cex = fscal, col=c("black", "gray40"))
    dev.off();
    
}

