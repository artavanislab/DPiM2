#!/pkg/R/3.1.2/EL6/x86_64/bin/Rscript

##############################################################################80

{
    ##add nrBaitWinner to the apms data
    ## then break the file up into digestible chunnks

    annotation <- read.table("/home/glocke/DPiM/augRemap/apmsData/annotation/metaDataTable.tap.geneSymbol_RAO_GLL_112116.txt", quote="", sep="\t", header = T)
    apms <- read.table("/home/glocke/DPiM/augRemap/apmsData/augRemap.dp4.newBait.rmDupes.pepFDR.sumIso.applyLC.untrans.tscCut.rebrand", header = T)
    ann2 <- annotation[,c("search_id", "nrBaitWinner")]
    apms2 <- merge(apms, ann2, by.x="search_id", by.y="search_id")
    
    expts <- levels(as.factor(apms2$search_id))
    groupBy <- 1000
    while (i < length(expts)) {
    rwz <- (i-1)*groupBy
}


{
    library(igraph)
    netTab <- read.table("/home/glocke/DPiM/augRemap/nrBait_yRej_nBrand_11-16-2016/nrBait.net", header = T)
    net <- graph_from_data_frame(netTab, directed=F)
    comm <- lapply(c(2:10), function(s) { cluster_walktrap(net, steps=s) })

    {
        ## takes 4 eva!
        connComp <- decompose(net)
        vcount(connComp[1])
        comm2 <- cluster_spinglass(connComp[[1]], spins=25)
    }

    {
        ## takes 4 eva!
        ##connComp <- decompose(net)
        ##vcount(connComp[1])
        comm2 <- cluster_leading_eigen(connComp[[1]])
        ##comm2 <- cluster_leading_eigen(net)
    }
    {
        ##connComp <- decompose(net)
        ##vcount(connComp[1])
        comm2 <- cluster_louvain(connComp[[1]])
        ##comm2 <- cluster_leading_eigen(net)
    }
}

{
    net <- read.table("/home/glocke/DPiM/augRemap/nrBait_yRej_nBrand_11-16-2016/nrBait.net.i2.annotated.txt", header = T, sep="\t", quote="")
    net$anyShared <- ! (is.na(net$sharedReactome) & is.na(net$sharedGO) & is.na(net$sharedPfam))
    curv <- sapply(c(1:nrow(net)), function(i) {
                       1 - sum(net$anyShared[c(1:i)])/i
                   })
    plot(1-curv, type="l", log="x", xlim=c(10, nrow(net)))
    plot(1-curv, type="l", main="Fraction of edges sharing an annotation",
         xlab="Edge rank", ylab="Fraction among edges up to rank")
}

{
    ## find sid's not in the bait id table
    statFile <- '/home/glocke/DPiM/augRemap/apmsData/augRemap.dp4.statsBySID'
    baitFile <- '/home/glocke/DPiM/augRemap/apmsData/annotation/baitIdentityTable.11-07-2016.tsv'
    sidTab <- read.table(statFile, header = T)
    baitTab <- read.table(baitFile, header = T)
    subset(sidTab, ! search_id %in% baitTab$search_id)
}

{
    ## search for fbgn's not in aveLen
    library(org.Dm.eg.db)

    aveLenFile <- '/home/glocke/DPiM/nsaf/REVdmel-all-translation-r6.07_TAGS_sorted_trEl_vir.aveLen.tsv'
    baitFile <- '/home/glocke/DPiM/augRemap/apmsData/annotation/baitIdentityTable.11-08-2016.tsv'
    outFile <- '/home/glocke/DPiM/augRemap/apmsData/annotation/notFoundBaits.11-08-2016.tsv'
    bobFile <- '/home/glocke/DPiM/augRemap/apmsData/annotation/bobsUpdate.tsv'
    lenTab <- read.table(aveLenFile, header = T)
    baitTab <- read.table(baitFile, header = T)
    bobTab <- read.table(bobFile, header = T)
    
    kept <- subset(baitTab, retain == 'Retain')
    notFound <- subset(kept, ! bait %in% lenTab$fbgn)
    notFoundBait <- as.character(unique(notFound$bait))
    
    notFoundTab <- select(org.Dm.eg.db, notFoundBait,
                          c("FLYBASE", "ENTREZID", "SYMBOL"),
                          keytype="FLYBASE")

    write.table(notFoundTab, file=outFile, quote=F, sep="\t", row.name=F)

    subset(bobTab, ! SubmittedID %in% lenTab$fbgn)
}


multiGrep <- function(patterns, vect, exact=TRUE) {
    greppp <- function(p) { grep(p, vect) }
    if (exact) {
        greppp <- function(p) { which(vect == p) }
    }
    sapply(patterns, greppp)
}

{
    args <- commandArgs(TRUE)
    x <- as.numeric(args[1])
    print(x*10)
    
    
    zscore <- function(x, bg) {
        mn <- mean(bg)
        ss <- sd(bg)
        return((x-mn)/ss)
    }
}
