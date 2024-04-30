##012345678901234567890123456789012345678901234567890123456789012345678901234567

library(biomaRt)

{
    dmel <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                    dataset="dmelanogaster_gene_ensembl",
                    host = "jul2015.archive.ensembl.org")
    outDir <- '/home/glocke/DPiM/uniprot/'
}

{
    raw <- read.table("/home/glocke/DPiM/dpim4/withInstr/apmsData/all.combined.05-04-2016.newFBgn.rmDupes.sumIso", header = T)
    allGenes <- union(raw$prey_ref, raw$bait_ref)
    trans <- getBM(attributes = c("flybase_gene_id", "uniprot_sptrembl"),
                   filters = "flybase_gene_id", values = allGenes, mart = dmel)
    write.table(trans, file=paste(outDir, "/biomartTrans.txt", sep=""),
                quote=F, row.names=F, sep="\t")
    write(allGenes, file=paste(outDir, "/allProt.tsv", sep=""))
}

eids <- c("35851", "42737", "34411", "35851", "33196", "35851", "43767", "34085", "41565", "41965", "36382", "32772", "41565", "34411", "37618", "43767", "35851", "41565", "40701", "40701", "37618", "33196", "42737", "42737", "42737")

# filter entrezgene
# attributes entrezgene, flybase_gene_id

eid2fb <- getBM(attributes=c("flybase_gene_id"), filters=c("entrezgene"), values=eids, mart=dmel)
