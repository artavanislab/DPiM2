## return a data frame with columns:
## protID, TSC, numExp, numBait
## TSC is the sum of all total_peptides for this protein (where bait!=prey)
## numExp is the number of search_ids containing this protein
## numBait is the number of search_id's with this protein as a bait

safe <- apms
apms <- apms[1:1000,]
apms <- safe

ap2 <- apms[1:10,]
ll <- levels(ap2$prey_ref)

apms[apms$prey_ref == 1]
apms$prey_ref[1:1000]
subset(apms, prey_ref == "FBgn0263391")
table(apms$prey_ref)

class(apms$prey_ref)

broadProteinStats <- function(apms) {
    allPR <- unique(apms$prey_ref)
    ##allPID <- unique(c(allBR, allPR))

    
    # get TSC for every prey
    tscP <- aggregate(apms[,"total_peptides"], by=list(apms$prey_ref),
                        FUN=sum)
    row.names(tscP) <- tscP$Group.1
    baitbait <- subset(apms, as.character(bait_ref) == as.character(prey_ref))
    bb_agg <- aggregate(baitbait[,"total_peptides"], by=list(baitbait$prey_ref),
                        FUN=sum)
    row.names(bb_agg) <- bb_agg$Group.1
    for (b in bb_agg$Group.1) {
        tscP[b,]$x <- tscP[b,]$x - bb_agg[b,]$x
    }

    # get TSC for every bait
    tscB <- aggregate(apms[,"total_peptides"], by=list(apms$bait_ref),
                        FUN=sum)
    row.names(tscB) <- tscB$Group.1
    bb_agg <- aggregate(baitbait[,"total_peptides"], by=list(baitbait$prey_ref),
                        FUN=sum)
    row.names(bb_agg) <- bb_agg$Group.1
    for (b in bb_agg$Group.1) {
        tscB[b,]$x <- tscB[b,]$x - bb_agg[b,]$x
    }

    preyTab <- table(apms$prey_ref)
    baitTab <- table(factor(apms$bait_ref, levels=levels(apms$prey_ref)))

    ret <- as.data.frame(t(
        vapply(allPR, function(x) {
                   tscb <- tscB[x,]$x
                   if (is.na(tscb)) {
                       tscb <- 0
                   }
                   tscp <- tscP[x,]$x
                   if (is.na(tscp)) {
                       tscp <- 0
                   }
                   return(c(as.character(x), tscb, tscp, preyTab[x],
                            baitTab[x]))
               }, c("FBgn0011715", 1, 1, 1, 1))
        ))
    colnames(ret) <- c("protID", "baitTSC", "preyTSC", "preyBatches", "baitBatches")
    head(ret)

    return(ret)
}



apms <- read.table("~/DPiM/data/dpim3.12232014.nrBait",
                   col.names = c("tap_id", "search_id", "sample_date",
                       "total_peptides", "bait_ref", "prey_ref"))

apms$sample_date <- as.Date(apms$sample_date)

write.table(bps, file="globalProteinStats.nrBait.tab", row.names="F")
badTime <- subset(apms, sample_date > as.Date("2009-06-01") &
                      sample_date < as.Date("2010-06-01"))

bps <- broadProteinStats(apms)
bps_bad <- broadProteinStats(badTime)

head(bps)
head(bps_bad)

ap2 <- apms[1:1000, ]
bps <- broadProteinStats(ap2)
head(bps)

allBR <- as.character(unique(ap2$bait_ref))
allPR <- as.character(unique(ap2$prey_ref))
unique(c(allBR, allPR))

length(unique(ap2$prey_ref))
unique(ap2$prey_ref)

as.character(allPR)

unique(list(c(1:10), c(2:11)))

s <- subset(ap2, bait_ref == "I love peanuts")
length(unique(s$sample_id))
