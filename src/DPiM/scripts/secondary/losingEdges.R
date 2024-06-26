## seek nodes that with high degree in dpim2 that lose edges in dpim3

{
    library(igraph)
    dpim2_df <- read.table("/home/glocke/DPiM/data/dpim2_nrtap.120123")
    dpim3_df <- read.table("/home/glocke/DPiM/data/dpim3.12232014.nrBait.network", header = T)
    dpim2_df <- dpim2_df[c(5,6,1,2,3,4),]
    dpim2 <- graph_from_data_frame(dpim2_df, directed=F)
}
