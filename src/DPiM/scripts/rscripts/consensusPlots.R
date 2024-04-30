##############################################################################80
{
    cons1 <- read.table("/home/glocke/DPiM/dpim4/witInstr/cons1_01-30-2016.tsv",
                        header = T);
}


{

    tab <- sapply(table(cons1$support), function(x) return(x));
    pdf(file="~/DPiM/plots/cons1_01-30-2016_support.pdf", height=6, width=6)
    plot(tab, log="y", type="h", main="Histogram of edge support",
         xlab="Support", ylab="Count")
    dev.off()

    
    pdf(file="~/DPiM/plots/cons1_01-30-2016_supportECDF.pdf", height=6, width=6)
    plot(ecdf(cons1$support), main="Empirical CDF of edge support",
         xlab="Support", ylab="Percentile", yaxt="n")
    axis(2, at=c(0:5)/5, labels=c(0:5)*20, las=2)
    grid()
    dev.off();

}

{
    fc <- read.table("/home/glocke/DPiM/dpim4/witInstr/consensus1_01-30-2016/consensus1_01-30-2016.stats.foldChange.tsv", header = T)
    pdf(file="~/DPiM/plots/cons1_01-30-2016_foldChange.pdf", height=6, width=6)
    plot(fc$step, sqrt(fc$diff)/fc$score, xlim=c(0,500),
         xlab="Index", ylab="Fold change", main="Maximum fold change in HGScore per step")
    dev.off();
}

plot(fc$step, fc$step*sqrt(fc$diff)/(fc$score), xlim=c(0,500),
     xlab="Index", ylab="Fold change", main="Maximum fold change in HGScore per step")
