#!/usr/bin/env Rscript

##############################################################################80

## Take a file with one cluster per line (tab separated ids)

source("~/DPiM/scripts/rscripts/homebrew.R")


   ####### +  +  +  +  + ### 
  #####  Subroutines  #####  
 ### +  +  +  +  + #######   

collector <- function(egTab, goTab, reacTab, pfamTab) {
    reacTab <- merge(reacTab, egTab, by.x="ENTREZID", by.y="ENTREZID")

    ret <- rbindlist(list(goTab, reacTab[,c(idCol, "ENTREZID", "REID")],
                          pfamTab))
    ret <- as.data.frame(ret)
    names(ret) <- c(idCol, "ENTREZID", "ID")
    ret <- subset(ret, ! is.na(ret$ID))
    ret <- subset(ret, ret$ID != "RE:NA")

    ret <- ret[order(ret[,idCol]),]
    
    return(ret)
}

pfamTable <- function(allEG, pfamFile) {
    pfam <- read.table(pfamFile)
    colnames(pfam) <- c('uniprot', 'alignment_start', 'alignment_end', 'envelope_start', 'envelope_end', 'hmm_acc', 'hmm_name', 'type', 'hmm_start', 'hmm_end', 'hmm_length', 'bit_score', 'E-value', 'clan')
    pfam2 <- pfam[,qw("uniprot hmm_acc")]

    colz <- unique(c(idCol, qw("ENTREZID UNIPROT")))
    uniprotTab <- select(eg.db, keys=allEG[,idCol], columns=colz, keytype=idCol)
    pfam3 <- subset(pfam2, uniprot %in% uniprotTab$UNIPROT)

    ret <- merge(uniprotTab, pfam3, by.x="UNIPROT", by.y="uniprot")
    return(unique(ret[,qw("FLYBASE ENTREZID hmm_acc")]))
}

## EASE algorithm: hypergeometric test where you subtract one from the count
## args:
##   found = number of annotated nodes found
##   nAnnotated = number of annotated nodes in universe
##   nAll = number of nodes in universe
##   pulls = number of nodes in set to be annotated
EASE <- function(found, nAnnotated, nAll, pulls) {
    test <- found - 1 ## EASE's correction to fisher/HG test
    white <- nAnnotated
    black <- nAll - nAnnotated
    return(phyper(test-1, white, black, pulls, lower.tail=F))
}

hypGeoTest <- function(found, nAnnotated, nAll, pulls) {
    test <- found ## skip the correction
    white <- nAnnotated
    black <- nAll - nAnnotated
    return(phyper(test-1, white, black, pulls, lower.tail=F))
}

jaccard <- function(l1, l2) {
    inter <- length(intersect(l1, l2))
    union <- length(union(l1, l2))
    return(inter/union)
}

jaccardCollapse <- function(byTerm, jaccardMax = 0.9) {
    membrs <- lapply(byTerm, function(df) {
                         df <- data.frame(df)
                         return(df[,idCol])
                      })

    jacc <- function(namePair) { jaccard(membrs[[namePair[1]]], membrs[[namePair[2]]]) }
    print("mapping jaccard distance between gene sets, this should take about 15 minutes")
    memPairs <- combn(names(membrs), 2, simplify=F)
    jaccDist <- lapply(memPairs, jacc)


    ##vJacc <- Vectorize(jacc)
    ##system.time(memPairs <- combn(names(membrs), 2, simplify=F))
    ##   user  system elapsed 
    ## 10.034   0.000  10.037 
    ##system.time(jaccDist <- vJacc(memPairs))
    ##    user  system elapsed 
    ## 872.254   0.000 872.443 
    ##system.time(jaccDist <- lapply(memPairs, jacc))
    ##    user  system elapsed 
    ## 829.557   0.015 829.893
    ##
    ## only 14 minutes for 5 million pairs; not bad!
}

oneClustTest <- function(cl, goodTerms, nTerms, termNames) {
    ## for one cluster, find the most significant enrichment
    thisDF <- subset(goodTerms, ENTREZID %in% cl)
    thisTab <- table(thisDF$ID)
    thisTab <- thisTab[thisTab>1]
    if (length(thisTab) == 0 || nrow(thisDF) == 0) {
        return(data.frame(minP = 1, minQ = 1, term = "NA", sig=FALSE,
                          size = length(cl), bestCnt = 0, allCnt=NA, nWayTie=NA,
                          name=NA))
    }
    cmpTab <- termTab[names(thisTab)]
    ##pvals <- EASE(thisTab, cmpTab, nNodes, length(cl))
    pvals <- hypGeoTest(thisTab, cmpTab, nNodes, length(cl))
    qvals <- p.adjust(pvals, method="holm", n=nTerms)
    x <- which(qvals == min(qvals))
    nBest <- length(x) ## how many 
    if(1 < nBest) {
        x <- x[1]
    }
    bestTerm <- names(x)
    return(data.frame(minP = pvals[x], minQ = qvals[x], term = names(x),
                      sig=qvals[x]<0.05, size = length(cl), 
                      bestCnt = thisTab[x], allCnt = cmpTab[x],
                      nWayTie = nBest, name=termNames[bestTerm,"name"]))
}

   ####### +  +  +  +  + ### 
  #####   Main Code   #####  
 ### +  +  +  +  + #######   

clustFile <- ""
outFile <- ""
modes <- c("fly", "human")

## for now, these will have to remain hard-coded file names
defaults <- list(
    pfamFile = '/home/glocke/DPiM/pfam/uniprot2pfam_Dmelanogaster_7227.tsv',
    pfamNameFile = '/home/glocke/DPiM/pfam/Pfam-A.hmm.dat.my.tsv'
    )

standardExcludedTerms <- qw('GO:0003674 GO:0005515 GO:0005575 GO:0005634
GO:0005737 GO:0005739 GO:0005829 GO:0008150 GO:0005829 GO:0008150')
excludedTerms <- c()

{
    ## collect command line arguments
    argv <- commandArgs(TRUE)
    if (length(argv) < 2) {
        tmpModes <- modes
        tmpModes[1] <- paste("*", tmpModes[1], "*", sep="")
        modeString <- paste(tmpModes, collapse=",")
        write(paste("usage: clustPerLine.txt annotation.out < ",
                    modeString, " exclusions >", sep=""), stderr())
        stop("required args missing")
    }
    clustFile <- argv[1]
    outFile <- argv[2]

    if (length(argv) >= 3) {
        moode <- argv[3]
        if (! moode %in% modes) {
            stop(paste("don't understand the moode argument '", moode, "'",
                       sep=""))
        }
    } else {
        moode <- 'fly'
    }
    if (length(argv) >= 4) {
        specialVal = 'exclusions'
        if (argv[4] == specialVal) {
            excludedTerms <- standardExcludedTerms
        } else {
            write(paste("found argv[4] == '", argv[4], "', so doing nothing. ",
                        "for standard exclusions, set argv[4] to '", specialVal,
                        "'", sep=''), stderr())
        }
    }
}
minCluster <- 3 ## only compute enriched percentage for clusters of this size or greater

if (moode == 'human') {
    stop("you'll have some work to do to check into what's happening with idCol, ENTREZID")
}


library(org.Dm.eg.db)
library(org.Hs.eg.db)
library(GO.db)
library(reactome.db)
##library(PFAM.db)
library(data.table)

if (moode == 'fly') {
    eg.db <- org.Dm.eg.db
    idCol <- "FLYBASE"
} else if (moode == 'human') {
    eg.db <- org.Hs.eg.db
    idCol <- "ENTREZID"
} else {
    stop(paste("unknown moode '", moode, "'", sep=""))
}

{
    clusters <- readLines(clustFile)
    clusters <- lapply(clusters, qw)
    allIDs <- unlist(clusters)
    nNodes <- length(allIDs)
    allEG <- select(eg.db, allIDs, c(idCol, "ENTREZID"), keytype=idCol)
    allEG <- dropRowsWithNaN(allEG, quant=FALSE)
    row.names(allEG) <- allEG[,idCol]
    
    nativeCl <- clusters
    egCl <- lapply(clusters, function(cl) { ## cluster members in entrez
                       return(allEG[cl,"ENTREZID"])
                   })
}

{
    ## annotate the universe!  WOWWIE!
    allGO <- select(eg.db, allIDs, c(idCol, "ENTREZID", "GO"), keytype=idCol)
    allGO <- unique(allGO[,c(idCol, "ENTREZID", "GO")])
    allGO <- subset(allGO, ! GO %in% excludedTerms)
    allGO <- subset(allGO, ! is.na(GO))
    termNames <- Term(GOTERM)
    goNames <- data.frame(term = names(termNames), name=termNames)

    allReac <- select(reactome.db, allEG$ENTREZID, c("ENTREZID", "PATHID"))
    allReac <- subset(allReac, ! is.na(PATHID))
    reacNames <- select(reactome.db, unique(allReac$PATHID),
                        c("PATHID", "PATHNAME"), keytype="PATHID")
    reacNames$MYNAME <- gsub("Drosophila melanogaster: ", "",
                             reacNames$PATHNAME)

    allReac$REID <- sapply(allReac$PATHID, function(id) {
                               paste("RE:", id, sep="")
                           })
    reacNames$REID <- sapply(reacNames$PATHID, function(id) {
                                 paste("RE:", id, sep="")
                             })
    row.names(reacNames) <- reacNames$REID

    allPfam <- pfamTable(allEG, defaults$pfamFile);
    pfamNames <- read.table(defaults$pfamNameFile, header = T)

    allNames <- rbindlist(list(goNames,reacNames[,c("REID", "MYNAME")],
                               pfamNames[,qw("term uniqName")]))
    allNames <- as.data.frame(allNames)
    allNames$term <- as.character(allNames$term)
    allNames$name <- as.character(allNames$name)
    row.names(allNames) <- allNames$term

    ## collect all databases into a single data frame
    allTerms <- collector(allEG, allGO, allReac, allPfam)
    
    allTerms$name <- allNames[allTerms$ID, "name"];
    
    ## this would probably be more efficiently done using table than split
    listByTerm <- split(allTerms, allTerms$ID)
    sizes <- sapply(listByTerm, nrow)
    listByTerm <- listByTerm[sizes > 1]

    nTerms <- length(listByTerm)

    goodTerms <- subset(allTerms, ID %in% names(listByTerm))
    termTab <- table(goodTerms$ID)
    ##listByNode <- split(allTerms, allTerms$ENTREZID)
}

{
    testList <- lapply(egCl, oneClustTest, goodTerms=goodTerms, nTerms=nTerms, termNames=allNames)
    tests <- rbindlist(testList)
    tests$i <- 1:nrow(tests)

    min3Tests <- subset(tests, size>=3)
    enrichedPercent <- sprintf("%.1f", 100*sum(min3Tests$sig)/nrow(min3Tests))
    nEnriched <- sum(min3Tests$sig)
    
    write(paste("# tested", clustFile, "using GO.db and Reactome.db"),
          file=outFile)
    write(paste("# ", enrichedPercent, "% of clusters beat alpha = 0.05",
                sep=""), file=outFile, append=T)
    write(paste("#", nEnriched, "clusters beat alpha = 0.05"), file=outFile,
          append=T)
    write(paste("# performed p-value corrections for", nTerms, "tests"),
          file=outFile, append=T)
    write.table(tests, file=outFile, sep="\t", row.names=F, append=T)
    print(enrichedPercent)
}

