enList <- read.table("/home/glocke/DPiM/dpim3.1/mcode1_10-12-2015/all.enrich.list")
flList <- read.table("/home/glocke/DPiM/dpim3.1/mcode1_10-12-2015/all.flagNet.list")

pears <- sapply(1:nrow(enList), function(i) {
                    en <- read.table(enList$V1[i], header = T);
                    fl <- read.table(flList$V1[i], header = T);
                    return(cor(en$p, fl$p))
                }
                )

for (i in 1:nrow(enList)) {
    
}
