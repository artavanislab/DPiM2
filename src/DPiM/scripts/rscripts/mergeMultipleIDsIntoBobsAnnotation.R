#!/pkg/R/3.1.2/EL6/x86_64/bin/Rscript

##############################################################################80


source("~/DPiM/scripts/rscripts/homebrew.R")

multiGrep <- function(patterns, vect, exact=TRUE) {
    greppp <- function(p) { grep(p, vect) }
    if (exact) {
        greppp <- function(p) { which(vect == p) }
    }
    sapply(patterns, greppp)
}

{
    multFile <- '/home/glocke/DPiM/augRemap/apmsData/annotation/metaDataTable3.tap.correct1.multi_RAO.txt'
    baseFile <- '/home/glocke/DPiM/augRemap/apmsData/annotation/metaDataTable.tap.geneSymbol_RAO_GLL_103116.txt'
    id2idFile <- '/home/glocke/DPiM/augRemap/apmsData/annotation/idString2SearchID.tsv'
    
    base <- read.delim(baseFile, header = T, sep="\t", stringsAsFactors=F)
    ## remove sample_date column, which excel has helpfully melted into slag
    {
        dateCol <- which(colnames(base)=="sample_date")
        base <- base[,-dateCol]
    }
    mult <- read.delim(multFile, header = T, sep="\t", stringsAsFactors=F)
    mult <- subset(mult, !is.na(search_id))

    id2id <- read.table(id2idFile, header = T, stringsAsFactors=F)
    
    ##> colnames(mult)
    ##[1] "search_id"      "ms_inst_run_id" "prevbait"       "tap_id"        
    ##[5] "newbait"        "sample_date"    "id_string"      "rejected"      
    ##[9] "wonky"          "tap_test"       "count"          "order_per_date"
    ##[13] "Comments"
    colnames(mult)[13] <- "multiComments"
    colnames(mult)[5] <- "Final_Updated_Bait_ID"
    ## amplify on the multiComments
    ordinaryMults <- which(nchar(mult$multiComments)==0)
    mult$multiComments[ordinaryMults] <- "multiple"

    dupeCols <- c( "Final_Updated_Bait_ID", "tap_id", "ms_inst_run_id", "rejected", "wonky")
    mult2 <- mult[,c("id_string", "order_per_date", "multiComments", dupeCols)]
    
    combn <- merge(base, mult2, by.x="id_string", by.y="id_string", all.x=T, all.y=T)##, incomparables=NA)
    combn <- merge(combn, id2id, by.x="id_string", by.y="id_string", all.x=T)##, incomparables=NA)
    nrow(combn)
    multiRows <- which(!is.na(combn$multiComments))

    ## check that there is a a search_id for each row
    noSID <- which(is.na(combn$search_id))
    if (length(noSID)) {
        errMsg <- paste("I found the following rows lacking a search_id: c(", paste(disagr, collapse=", "), ")")
        stop(errMsg)
    }
    
    ## check that new and old id's agree
    disagr <- which( (!is.na(combn$Final_Updated_Bait_ID.x)) &&
                        (!is.na(combn$Final_Updated_Bait_ID.y)) &&
                            combn$Final_Updated_Bait_ID.x !=
                                combn$Final_Updated_Bait_ID.x )
    if (length(disagr)) {
        errMsg <- paste("I found the following rows with different 'final' bait id's: c(", paste(disagr, collapse=", "), ")")
        stop(errMsg)
    }

}

{
    ##
    ## copy the information in duplicated columns from .y to .x where necessary
    ##

    ## first check that there are no disagreements
    for (d in dupeCols) {
        d.x <- paste(d, ".x", sep="")
        d.y <- paste(d, ".y", sep="")
        disagr <- which(!is.na(combn[multiRows, d.x]) &&
                            combn[multiRows, d.x] != combn[multiRows, d.y])
        if (length(disagr)) {
            errMsg <- paste("I found the following rows where base and mult disagree on ", d, ": c(", paste(disagr, collapse=", "), ")")
            stop(errMsg)
        }
    }

    ## now copy 
    dupeCols.x <- paste(dupeCols, ".x", sep="")
    dupeCols.y <- paste(dupeCols, ".y", sep="")

    combn[multiRows, dupeCols.x] <- combn[multiRows, dupeCols.y]

    ## remove one of the duplicate columns and rename the other
    dupeColsN.x <- multiGrep(dupeCols.x, colnames(combn))
    dupeColsN.y <- multiGrep(dupeCols.y, colnames(combn))
    colnames(combn)[dupeColsN.x] <- dupeCols
    combn <- combn[,-dupeColsN.y]
    
}

{
    ## say the retain information
    combn$RAO_Determination[multiRows] <- "Retain"
    ## a few of these multiple rids were actually duplicates, i.e. multiple
    ## rows in the table corresponding the same occurrence in meat-space
    badMults <- grep("This ID does not correspond", combn$multiComments)
    combn$RAO_Determination[badMults] <- "Remove"
    ## RAO: I have all 283 runs from May 14, 2014 listed as "Failed Runs.
    badDay <- grep("2014-05-14", combn$multiComments)
    combn$RAO_Determination[badDay] <- "Remove"
    ## for debuging purposes
    combn$RAO_Confidence[multiRows] <- "High"
}

{
    wholeTableFile <- "/home/glocke/DPiM/augRemap/apmsData/annotation/metaDataTable.tap.geneSymbol_RAO_GLL_110716.txt"
    write.table(combn, file=wholeTableFile, sep="\t", row.names=F, quote=F)

    baitKeyFile <- "/home/glocke/DPiM/augRemap/apmsData/annotation/baitIdentityTable.11-07-2016.tsv"
    baitKey <- combn[,c("search_id", "Final_Updated_Bait_ID",
                        "RAO_Determination")]
    colnames(baitKey) <- qw("search_id       bait    retain")
    write.table(baitKey, file=baitKeyFile, sep="\t", row.names=F, quote=F)
}


if (0) {
    ## some updates were made by hand
    wholeTableFile <- "/home/glocke/DPiM/augRemap/apmsData/annotation/metaDataTable.tap.geneSymbol_RAO_GLL_110716.txt"
    combn <- read.delim(wholeTableFile, header = T, sep="\t")

    baitKeyFile <- "/home/glocke/DPiM/augRemap/apmsData/annotation/baitIdentityTable.11-07-2016.tsv"
    baitKey <- combn[,c("search_id", "Final_Updated_Bait_ID", "RAO_Determination")]
                        
    colnames(baitKey) <- qw("search_id       bait    retain")
    write.table(baitKey, file=baitKeyFile, sep="\t", row.names=F, quote=F)
}
