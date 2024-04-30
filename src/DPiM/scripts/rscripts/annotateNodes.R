#!/usr/bin/env Rscript

##############################################################################80

# KJ: this is a direct copy of the annotateClusters.R file
# except only write out the allTerms variable content



## Take a file with one cluster per line (tab separated ids)

source("~/DPiM/scripts/rscripts/homebrew.R")


   ####### +  +  +  +  + ### 
  #####  Subroutines  #####  
 ### +  +  +  +  + #######   

#collector <- function(egTab, goTab, reacTab, pfamTab) {
#    reacTab <- merge(reacTab, egTab, by.x="ENTREZID", by.y="ENTREZID")
#
#    ret <- rbindlist(list(goTab, reacTab[,c(idCol, "ENTREZID", "REID")],
#                          pfamTab))
#    ret <- as.data.frame(ret)
#    names(ret) <- c(idCol, "ENTREZID", "ID")
#    ret <- subset(ret, ! is.na(ret$ID))
#    ret <- subset(ret, ret$ID != "RE:NA")
#
#    ret <- ret[order(ret[,idCol]),]
#    
#    return(ret)
#}

collector <- function(egTab, goTab, reacTab, pfamTab, fbTab, otTab, mbTab) {
    reacTab <- merge(reacTab, egTab, by.x="ENTREZID", by.y="ENTREZID")

    ret <- rbindlist(list(goTab, reacTab[,c(idCol, "ENTREZID", "REID")],
                          pfamTab, fbTab, otTab, mbTab))
    ret <- as.data.frame(ret)
    names(ret) <- c(idCol, "ENTREZID", "ID")
    ret <- subset(ret, ! is.na(ret$ID))
    ret <- subset(ret, ret$ID != "RE:NA")

    ret <- ret[order(ret[,idCol]),]

    return(ret)
}



pfamTable <- function(allEG, pfamFile) {
    #pfam <- read.table(pfamFile)
    #colnames(pfam) <- c('uniprot', 'alignment_start', 'alignment_end', 'envelope_start', 'envelope_end', 'hmm_acc', 'hmm_name', 'type', 'hmm_start', 'hmm_end', 'hmm_length', 'bit_score', 'E-value', 'clan')
    #pfam2 <- pfam[,qw("uniprot hmm_acc")]

	pfam <- read.table(pfamFile, header=T, sep="\t", strip.white=T, stringsAsFactors=F, quote="")
	pfam2 <- pfam[,qw("UNPACCESSION PFAM_ACC")]

    colz <- unique(c(idCol, qw("ENTREZID UNIPROT")))
    uniprotTab <- select(eg.db, keys=allEG[,idCol], columns=colz, keytype=idCol)
    #pfam3 <- subset(pfam2, uniprot %in% uniprotTab$UNIPROT)
	pfam3 <- subset(pfam2, UNPACCESSION %in% uniprotTab$UNIPROT)

    #ret <- merge(uniprotTab, pfam3, by.x="UNIPROT", by.y="uniprot")
	ret <- merge(uniprotTab, pfam3, by.x="UNIPROT", by.y="UNPACCESSION")
    #return(unique(ret[,qw("FLYBASE ENTREZID hmm_acc")]))
	return(unique(ret[,qw("FLYBASE ENTREZID PFAM_ACC")]))
}

# Open Targets genetic association table
# /home/kli3/proj/OpenTargets/symbol2diseaseAssociation.opentargets.csv
otTable <- function(allEG, otFile) {

    ot <- read.table(otFile, header=T, sep=",", strip.white=T, quote="\"", stringsAsFactors=F)
	# split the disease terms
	s <- strsplit(ot$n.disease_opentarget, split = "\\|\\|")
	ot2 <- data.frame(GENE = rep(ot$n.symbol, sapply(s, length)), DISEASE = unlist(s))
	ot2 <- unique(ot2)
    ot3 <- data.frame(OpenTargets_ID="", term=unique(ot2$DISEASE))
	ot3$OpenTargets_ID <- paste("OpenTargets_ID:",rownames(ot3),sep="")
    ot3$uniqName <- paste(ot3$term, ot3$OpenTargets_ID, sep=" ")

	ot4 <- merge(ot2, ot3, by.x="DISEASE", by.y="term")

	# diopt as flybase to human translation /home/glocke/DPiM/corum/DIOPT_Translations_for_George_110916.txt
	# TODO: check to see which column is the right FBgn id "FlyBaseID" or "FBgn.Search.Term", using "FlyBaseID" now
	# 05082017 we are going to "FlyBaseID" and ONLY use 'high' or 'moderate' confidence orthologs 
	diopt <- read.table("/home/glocke/DPiM/corum/DIOPT_Translations_for_George_110916_removed_97_special_char.txt", sep="\t", header=T, comment.char="#", quote="\"", stringsAsFactors=F)	
	
	# 072019 team discussed to include FBgns with only one human ortholog even it is below our cutoff
	# 072019 KJ added the code
	tres <- sort(table(diopt$FBgn.Search.Term))
	ones <- diopt[ which(diopt$FBgn.Search.Term %in% names(tres[tres==1])),] # FBgns only has 1 human ortholog
	ones <- ones[ones$Rank=='low'|ones$Rank=='None',] # rescuing low and None ranked ones

	diopt <- subset(diopt, Rank=='high'|Rank=='moderate',)
	
	# 072019 KJ added the code
	diopt <- rbind(diopt, ones) 	

	colz <- unique(c(idCol, qw("ENTREZID")))
    otTab <- select(eg.db, keys=allEG[,idCol], columns=colz, keytype=idCol)

	# first lookup FBgn in the diopt table
    ret <- merge(otTab, diopt, by.x="FLYBASE", by.y="FlyBaseID")
	# then use Human symbol to look up open targets table
	ret <- merge(ret, ot4, by.x="Human.Symbol", by.y="GENE")


    return(list(unique(ret[,qw("FLYBASE ENTREZID OpenTargets_ID")]), ot3))
}

# MetaBase disease association table
# /home/kli3/proj/TR/metabase/disease_associations.txt
mbTable <- function(allEG, mbFile) {
	
	mb <- read.table(mbFile, header=T, sep="\t", strip.white=T, quote="\"", stringsAsFactors=F)
    # split the disease terms

    mb$mb_ID <- paste("MetabaseID:",mb$disease_id,sep="")
	mb$uniqName <- paste(mb$disease_name, mb$mb_ID, sep=" ")

    mb <- subset(mb, !is.na(genesymbol), ) # remove rows with genesymbol as NA

    # diopt as flybase to human translation /home/glocke/DPiM/corum/DIOPT_Translations_for_George_110916.txt
    # 05082017 we are going to "FlyBaseID" and ONLY use 'high' or 'moderate' confidence orthologs
    diopt <- read.table("/home/glocke/DPiM/corum/DIOPT_Translations_for_George_110916_removed_97_special_char.txt", sep="\t", header=T, comment.char="#", quote="\"", stringsAsFactors=F)

    # 072019 team discussed to include FBgns with only one human ortholog even it is below our cutoff
    # 072019 KJ added the code
    tres <- sort(table(diopt$FBgn.Search.Term))
    ones <- diopt[ which(diopt$FBgn.Search.Term %in% names(tres[tres==1])),] # FBgns only has 1 human ortholog
    ones <- ones[ones$Rank=='low'|ones$Rank=='None',] # rescuing low and None ranked ones

    diopt <- subset(diopt, Rank=='high'|Rank=='moderate',)

    # 072019 KJ added the code
    diopt <- rbind(diopt, ones)

    colz <- unique(c(idCol, qw("ENTREZID")))
    mbTab <- select(eg.db, keys=allEG[,idCol], columns=colz, keytype=idCol)

    # first lookup FBgn in the diopt table
    ret <- merge(mbTab, diopt, by.x="FLYBASE", by.y="FlyBaseID")
    # then use Human symbol to look up metabase table
    ret <- merge(ret, mb, by.x="Human.Symbol", by.y="genesymbol")

    return(list(unique(ret[,qw("FLYBASE ENTREZID mb_ID")]), unique(mb[,qw("mb_ID uniqName")])))
}

# FlyBase disease association table (DOID)
# /home/glocke/DPiM/augRemap/apmsData/annotation/allele_human_disease_model_data_fb_2017_02.txt
fbTable <- function(allEG, fbFile) {
    doid <- read.table(fbFile, header=T, sep="\t", comment.char="#", quote="")
    doid2 <- doid[,qw("FBgn_ID DOID_ID")]
	
	colz <- unique(c(idCol, qw("ENTREZID")))
    fbTab <- select(eg.db, keys=allEG[,idCol], columns=colz, keytype=idCol)

	# TODO: question about the FBgn id version, check with Bob and George
	# 0508017 due to the small changes from version to version
	# we decided to ignore the difference for now, if needed, we will revisit this
	ret <- merge(fbTab, doid2, by.x="FLYBASE", by.y="FBgn_ID")

    return(unique(ret[,qw("FLYBASE ENTREZID DOID_ID")]))
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

#KJ modified the following function to collect top terms to cluster name
oneClustTest <- function(cl, goodTerms, nTerms, termNames) {
    #cat(paste(cl,collapse=":::"))
    #cat("\n")
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
    #x <- which(qvals == min(qvals))
    #nBest <- length(x) ## how many
    #if(1 < nBest) {
    #    x <- x[1]
    #}
    x <- sort(qvals, decreasing=F, index.return=T)
    bestTerm <- names(x$x[1:min(3,length(qvals))])
    return(data.frame(minP = pvals[x$ix[1]], minQ = qvals[x$ix[1]], term = paste(bestTerm, collapse=" ::: "),
                      sig=qvals[x$ix[1]]<0.05, size = length(cl),
                      bestCnt = thisTab[x$ix[1]], allCnt = cmpTab[x$ix[1]],
                      nWayTie = NA, name=paste(termNames[bestTerm,"name"], collapse=" ::: ")))
}

   ####### +  +  +  +  + ### 
  #####   Main Code   #####  
 ### +  +  +  +  + #######   

clustFile <- ""
outFile <- ""
modes <- c("fly", "human")

## for now, these will have to remain hard-coded file names
defaults <- list(
    #pfamFile = '/home/glocke/DPiM/pfam/uniprot2pfam_Dmelanogaster_7227.tsv',
	pfamFile = '/home/glocke/DPiM/augRemap/apmsData/annotation/drome-pfams.txt',
    #pfamNameFile = '/home/glocke/DPiM/pfam/Pfam-A.hmm.dat.my.tsv'
	pfamNameFile = '/home/glocke/DPiM/augRemap/apmsData/annotation/drome-pfams.txt',
	fbFile = '/home/glocke/DPiM/augRemap/apmsData/annotation/allele_human_disease_model_data_fb_2017_02.txt',
	otFile = '/home/kli3/proj/OpenTargets/symbol2diseaseAssociation.opentargets.csv',
	mbFile = '/home/kli3/proj/TR/metabase/disease_associations.txt'
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

#minCluster <- 3 ## only compute enriched percentage for clusters of this size or greater
#KJ made the change of minimum side of cluster to be considered as 2 per DPiM team meeting 04/28/2017
minCluster <- 2 ## only compute enriched percentage for clusters of this size or greater

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
    commLines <- grep("^#", clusters, perl=T) 
    clusters <- clusters[-commLines]
    clusters <- lapply(clusters, qw)
    clusterNames <- sapply(clusters, head, n=1)
    clusters <- lapply(clusters, function(cl) { cl[-1] })
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
    #pfamNames <- read.table(defaults$pfamNameFile, header = T)
	pfamNames <- read.table(defaults$pfamFile, header=T, sep="\t", strip.white=T, comment.char="#", stringsAsFactors=F, quote="")
	pfamNames$uniqName <- paste(pfamNames$DESCRIPTION,pfamNames$PFAM_ACC,sep=" ")


	## adding Open targets, flybase and metabase disease associations here
	# Flybase disease
	allFb <- fbTable(allEG, defaults$fbFile);
	fbNames <- read.table(defaults$fbFile, header=T, sep="\t", comment.char="#", quote="")
	fbNames$uniqName <- paste(fbNames$DOID_term, fbNames$DOID_ID, sep=" ")
	# Open Targets
	allOt <- otTable(allEG, defaults$otFile);
	otNames <- allOt[[2]]
	allOt <- allOt[[1]]
	#Metabase
	allMb <- mbTable(allEG, defaults$mbFile);
	mbNames <- allMb[[2]]
    allMb <- allMb[[1]]	

    #allNames <- rbindlist(list(goNames,reacNames[,c("REID", "MYNAME")],
    #                           pfamNames[,qw("term uniqName")]))
	allNames <- rbindlist(list(goNames,reacNames[,c("REID", "MYNAME")],
								unique(pfamNames[,qw("PFAM_ACC uniqName")]),
								unique(fbNames[,qw("DOID_ID uniqName")]),
								unique(otNames[,qw("OpenTargets_ID uniqName")]),
								unique(mbNames)
							  )
						 )
	
    allNames <- as.data.frame(allNames)
    allNames$term <- as.character(allNames$term)
    allNames$name <- as.character(allNames$name)
    row.names(allNames) <- allNames$term

    ## collect all databases into a single data frame
    #allTerms <- collector(allEG, allGO, allReac, allPfam)
	allTerms <- collector(allEG, allGO, allReac, allPfam, allFb, allOt, allMb)
    
	#KJ allTerms$name <- allNames[allTerms$ID, "name"];
    allTerms$name <- allNames[as.character(allTerms$ID), "name"];    
	write.table(allTerms, file=outFile, sep="\t", row.names=F, quote=F)
	quit("no") # exit without running the rest of the code

    ## this would probably be more efficiently done using table than split
    ## but i did it this way in anticipation of using jaccard to merge terms
    listByTerm <- split(allTerms, as.character(allTerms$ID))
    sizes <- sapply(listByTerm, nrow)
    listByTerm <- listByTerm[sizes > 1]

    nTerms <- length(listByTerm)

    goodTerms <- subset(allTerms, ID %in% names(listByTerm))
    termTab <- table(goodTerms$ID)
    ##listByNode <- split(allTerms, allTerms$ENTREZID)
}

{
    ## test all clusters
    testList <- lapply(egCl, oneClustTest, goodTerms=goodTerms, nTerms=nTerms, termNames=allNames)
    tests <- rbindlist(testList)
    tests <- cbind(clusterNames, tests)
    names(tests)[1] <- 'clustID'

#    min3Tests <- subset(tests, size>=3)
#    enrichedPercent <- sprintf("%.1f", 100*sum(min3Tests$sig)/nrow(min3Tests))
#    nEnriched <- sum(min3Tests$sig)
#	KJ modified the code to allow using minCluster
	minTests <- subset(tests, size>=minCluster)
    enrichedPercent <- sprintf("%.1f", 100*sum(minTests$sig)/nrow(minTests))
    nEnriched <- sum(minTests$sig)

}

{
    ## report
    tests <- tests[order(-tests$size),]
    tests <- format(tests, digits=5)

    write(paste("# tested", clustFile, "using GO.db, Reactome.db, the Pfam table provided by Krishna, Flybase disease association provide by Bob, OpenTargets disease association provided by Dan Park and Metabase disease asscoation pulled from Metabase 6.28"),
          file=outFile)
    write(paste("# ", enrichedPercent, "% of clusters beat alpha = 0.05",
                sep=""), file=outFile, append=T)
    write(paste("#", nEnriched, "clusters beat alpha = 0.05"), file=outFile,
          append=T)
    write(paste("# performed p-value corrections for", nTerms, "tests"),
          file=outFile, append=T)
    write.table(tests, file=outFile, sep="\t", row.names=F, append=T, quote=F)
    print(enrichedPercent)
}

