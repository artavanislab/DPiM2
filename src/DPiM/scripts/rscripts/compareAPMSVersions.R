expts <- read.table("/home/glocke/DPiM/newDpim2/compareVersions.apms.tsv.sym", header = T)

ids <- unique(expts$ms_inst_run_id )

noBait <- subset(expts, as.character(bait)!= as.character(prey))
{
    {
        cors <- sapply(ids, function(id) {
                           exp <- subset(noBait, ms_inst_run_id == id);
                           return(cor(exp$dp1TSC, exp$dp1aTSC))
                       })
        
        df <- data.frame(ms_inst_run_id = ids, cor=cors)
        df <- df[order(df$cor),]
        write.table(df, file="/home/glocke/DPiM/newDpim2/compareVersions.apms.cor_noBait.tab", quote=F, sep="\t", row.names=F)
    }

    {
        cors <- sapply(ids, function(id) {
                           exp <- subset(expts, ms_inst_run_id == id);
                           return(cor(exp$dp1TSC, exp$dp1aTSC))
                       })
        
        df <- data.frame(ms_inst_run_id = ids, cor=cors)
        df <- df[order(df$cor),]
        write.table(df, file="/home/glocke/DPiM/newDpim2/compareVersions.apms.cor.tab", quote=F, sep="\t", row.names=F)
    }
}

{
    corTab <- read.table("/home/glocke/DPiM/newDpim2/compareVersions.apms.cor.tab", header = T)
    ex2 <- merge(expts, corTab, by.x="ms_inst_run_id", by.y="ms_inst_run_id")
    ex2 <- ex2[order(ex2$cor),]
    ex2$cor <- round(ex2$cor, digits=4)

    ##corWithBait <- read.table("/home/glocke/DPiM/newDpim2/compareVersions.apms.cor_withBait.tab", header = T)
    ##c2 <- data.frame(ms_inst_run_id = corWithBait$ms_inst_run_id,
    ##corWithBait = corWithBait$cor)
    ##ex2 <- merge(ex2, c2, by.x="ms_inst_run_id", by.y="ms_inst_run_id")
    ##ex2$corWithBait <- round(ex2$corWithBait, digits=4)

    write.table(ex2, file="/home/glocke/DPiM/newDpim2/compareVersions.apms.txt", quote=F, sep="\t", row.names=F)
    
}

{
    sum(corTab$cor > 0.95)/nrow(corTab)
    sum(corTab$cor > 0.90)/nrow(corTab)
    sum(corTab$cor < 0.5)
    sum(corTab$cor < 0.5)/nrow(corTab)
    sum(corTab$cor > 0.5 & corTab$cor < 0.9)
    sum(corTab$cor > 0.5 & corTab$cor < 0.9)/nrow(corTab)
}
