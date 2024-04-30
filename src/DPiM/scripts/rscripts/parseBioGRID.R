{
    # biogrid tab2 format
    "
1       #BioGRID Interaction ID
2       Entrez Gene Interactor A
3       Entrez Gene Interactor B
4       BioGRID ID Interactor A
5       BioGRID ID Interactor B
6       Systematic Name Interactor A
7       Systematic Name Interactor B
8       Official Symbol Interactor A
9       Official Symbol Interactor B
10      Synonyms Interactor A
11      Synonyms Interactor B
12      Experimental System
13      Experimental System Type
14      Author
15      Pubmed ID
16      Organism Interactor A
17      Organism Interactor B
18      Throughput
19      Score
20      Modification
21      Phenotypes
22      Qualifications
23      Tags
24      Source Database
"
    # score has value '-' for all low throughput rows
}


lowT <- read.delim("/home/glocke/DPiM/biogrid/lowThru.tab", sep="\t")
fbgn <- read.table("/home/glocke/DPiM/flybase/fbgn_locusTag.tab", header = T)
fbgn <- subset(fbgn, locusTag != '-')
row.names(fbgn) <- fbgn$locusTag

lowT$fbgnA <- fbgn[lowT$Systematic.Name.Interactor.A, 'fbgn']
lowT$fbgnB <- fbgn[lowT$Systematic.Name.Interactor.B, 'fbgn']

lowT <- subset(lowT, lowT$Experimental.System.Type == 'physical')


edgeList <- data.frame(fbgnA = lowT$fbgnA, fbgnB = lowT$fbgnB, 
                       nameA = lowT$Official.Symbol.Interactor.A,
                       nameB = lowT$Official.Symbol.Interactor.B,
                       expSystem = lowT$Experimental.System,
                       author = lowT$Author,
                       pmid = lowT$Pubmed.ID, sourceDB = lowT$Source.Database,
                       bioGRIDID = lowT$X.BioGRID.Interaction.ID)

write.table(edgeList, file="/home/glocke/DPiM/biogrid/lowThruPhysical.tab",
            row.names=F)

library(igraph)

net <- graph_from_data_frame(edgeList)
conn <- clusters(
