argv <- commandArgs(TRUE)
if (length(argv) != 1) {
    stop("usage: clusterAnnotate.R <clusterPerLine.tsv> <parent.network>");
}
clustFile <- argv[1]
netFile <- argv[2]


library(org.Dm.eg.db)


clusterStrings <- scan(inFile, what=character())

clusters <- strsplit(clusterString, "[[:space:]]+")


