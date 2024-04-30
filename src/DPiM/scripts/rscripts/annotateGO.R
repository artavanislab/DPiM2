#!/usr/bin/env Rscript

##############################################################################80

## Take a file with one cluster per line (tab separated ids)
## Report only GO annotations evidenced by "Experimental Evidence codes"
## see http://www.geneontology.org/page/guide-go-evidence-codes

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

   ####### +  +  +  +  + ### 
  #####   Main Code   #####  
 ### +  +  +  +  + #######   

clustFile <- ""
outFile <- ""
modes <- c("fly", "human")

## for now, these will have to remain hard-coded file names
defaults <- list(
    pfamFile = '/home/glocke/DPiM/pfam/uniprot2pfam_Dmelanogaster_7227.tsv',
    pfamNameFile = '/home/glocke/DPiM/pfam/Pfam-A.hmm.dat.my.tsv',
    keep = 'EXP,IDA,IPI,IMP,IGI,IEP'
    )

{
    ## collect command line arguments
    argv <- commandArgs(TRUE)
    if (length(argv) < 2) {
        tmpModes <- modes
        tmpModes[1] <- paste("*", tmpModes[1], "*", sep="")
        modeString <- paste(tmpModes, collapse=",")
        write(paste("usage: clustPerLine.txt annotation.out < ",
                    modeString, " >", sep=""), stderr())
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
}

if (moode == 'human') {
    stop("you'll have some work to do to check into what's happening with idCol, ENTREZID")
}


library(org.Dm.eg.db)
library(org.Hs.eg.db)
library(GO.db)
##library(reactome.db)
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
    commLines <- grep("^#", clusters, perl=T) 
    clusters <- clusters[-commLines]
    clusters <- lapply(clusters, qw)
    clusterNames <- sapply(clusters, head, n=1)
    clusters <- lapply(clusters, function(cl) { cl[-1] })
    allIDs <- unlist(clusters)
    nNodes <- length(allIDs)
    allEG <- select(eg.db, allIDs, c(idCol, "ENTREZID", "SYMBOL"),
                    keytype=idCol)
    allEG <- dropRowsWithNaN(allEG, quant=FALSE)
    row.names(allEG) <- allEG[,idCol]

    keepCodes <- unlist(strsplit(defaults$keep, ','))
}

{
    ## annotate the universe!  WOWWIE!
    allGO <- select(eg.db, allIDs, c(idCol, "ENTREZID", "GO"), keytype=idCol)
    allGO <- subset(allGO, EVIDENCE %in% keepCodes)
    allGO <- unique(allGO[,c(idCol, "ENTREZID", "GO")])
    allGO <- subset(allGO, ! is.na(GO))
    termNames <- Term(GOTERM)
    goNames <- data.frame(term = names(termNames), name=termNames)

    
    allNames <- goNames
    allNames$term <- as.character(allNames$term)
    allNames$name <- as.character(allNames$name)
    row.names(allNames) <- allNames$term

    
    ## collect all databases into a single data frame
    allTerms <- allGO
    names(allTerms) <- c(idCol, "ENTREZID", "ID")    
    
    allTerms$name <- allNames[allTerms$ID, "name"];
    write.table(allTerms, file=outFile, sep="\t", row.names=F, quote=F)
}
