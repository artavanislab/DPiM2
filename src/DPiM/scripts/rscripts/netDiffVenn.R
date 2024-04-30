##345678901234567890123456789012345678901234567890123456789012345678901234567890

{
    library(VennDiagram)

    edgeList <- function(node1, node2) {
        df <- data.frame(n1=node1, n2=node2)
        apply(df, 1, function(row) {
                  paste(sort(row), collapse="_")
              })
    }
    edgeList2 <- function(df) {
        if (ncol(df) != 2) {
            stop("edgeList2 requires a data.frame with precisely two columns")
        }
        apply(df, 1, function(row) {
                  paste(sort(row), collapse="_")
              })
    }
}


{
    dpim1Tab <- read.table("/home/glocke/DPiM/prevDPIM/DPIM1_scores.r6.07.updateFBgn.tsv", header = T)
    dpim1JulTab <- read.table("/home/glocke/DPiM/newDpim1_2/newNrBait/newNrBait.net", header = T)
    dpim1Aug1Tab <- read.table("/home/glocke/DPiM/augRemap/nrBait11-08-2016/nrBait.net", header = T)
    edgeStrings <- list(
        dpim1=edgeList(dpim1Tab$Interactor_1, dpim1Tab$Interactor_2),
        dpim1Jul=edgeList(dpim1JulTab$protein1, dpim1JulTab$protein2),
        dpim1Aug1=edgeList(dpim1Aug1Tab$protein1, dpim1Aug1Tab$protein2)
        )
    
}

{
    vennFile = "/home/glocke/DPiM/plots/dpimVersions_11-08-2016.sharedEdges.tiff"
    mn <- list(dpim1="Published", dpim1Jul="DPIM1.2\nJuly 2016", dpim1Aug1="November 2016")
    venn.diagram(edgeStrings, filename=vennFile, sub.fontfamily="sans-serif",
                 main.fontfamily="sans-serif",
                 category.names = mn, fill=rainbow(3),
                 alpha=0.5)
}


{
    hgScore <- read.table("/home/glocke/DPiM/human/biopStyleData/nrBait_rebranded_09-26-2016/nrBait.net", header = T)
    comp85 <- read.table("/home/glocke/DPiM/human/CompPASS/BioPlex_8467map-edges.net", header = T)
    comp59 <- read.table("/home/glocke/DPiM/human/CompPASS/BioPlex.v5884_HCIs.net", header = T)
    frames <- list(c59=comp59[,c(1:2)], c85=comp85[,c(1:2)], hg=hgScore[,c(1:2)])
    edgeStrings <- lapply(frames, edgeList2);
}

{
    vennFile = "/home/glocke/DPiM/plots/hgVsCompPASS.sharedEdges_10-04-2016.tiff"
    mn <- list(c59="BioPlex 2.0", c85="BioPlex 2.5", hg="Human-HG (2.5 data)")
    venn.diagram(edgeStrings, filename=vennFile, sub.fontfamily="sans-serif",
                 main.fontfamily="sans-serif",
                 category.names = mn, fill=rainbow(3),
                 alpha=c(0.5, 0.5, 0.5))
}


{
    nodeStrings <- list(
        c59=unique(c(comp59$protein1, comp59$protein2)),
        c85=unique(c(comp85$protein1, comp85$protein2)),
        hg=unique(c(hgScore$protein1, hgScore$protein2))
        )
    vennFile = "/home/glocke/DPiM/plots/hgVsCompPASS.sharedNodes_10-04-2016.tiff"
    mn <- list(c59="BioPlex 2.0", c85="BioPlex 2.5", hg="Human-HG (2.5 data)")
    venn.diagram(nodeStrings, filename=vennFile, sub.fontfamily="sans-serif",
                 main.fontfamily="sans-serif",
                 category.names = mn, fill=rainbow(3),
                 alpha=c(0.5, 0.5, 0.5))
    
}

###
### old syntax

{
    minSTab <- read.table("/home/glocke/DPiM/dpim4/withInstr/consensus5_nets/support0001/support0001.net", header = T)
    maxSTab <- read.table("/home/glocke/DPiM/dpim4/withInstr/consensus5_nets/support2000/support2000.net", header = T)
    dpim1Tab <- read.table("/home/glocke/DPiM/prevDPIM/DPIM1_scores.updateFBgn.tsv", header = T)
    ##randTab <- read.table("/home/glocke/DPiM/dpim4/withInstr/consensus4_04-04-2016/nets/hgc4_001337.fdr.out")

    edgeStrings <- list(
        minS=edgeList(minSTab$protein1, minSTab$protein2),
        maxS=edgeList(maxSTab$protein1, maxSTab$protein2),
        dpim1=edgeList(dpim1Tab$Interactor_1, dpim1Tab$Interactor_2)
        )
    edgeStrings3 <- list(maxS = edgeStrings$maxS,
                         static=edgeList(staticTab$protein1, staticTab$protein2),
                         rand=edgeList(randTab$V1, randTab$V2))
}

{
    vennFile = "/home/glocke/DPiM/plots/consensus5_nets.sharedEdges.tiff"
    mn <- list(minS = "Inclusive", maxS="Exclusive", dpim1="DPIM1")
    venn.diagram(edgeStrings, filename=vennFile, sub.fontfamily="sans-serif",
                 main.fontfamily="sans-serif",
                 category.names = mn, fill=c("orangered", "magenta", "green"),
                 alpha=c(0.5, 0.5, 0.5))
}

{
    staticTab <- read.table("/home/glocke/DPiM/dpim4/withInstr/staticCore/static.net", header = T)
    support1800Tab <- read.table("/home/glocke/DPiM/dpim4/withInstr/consensus5_nets/support1800/support1800.net", header = T)
    edgeStrings2 <- list(cons = edgeList(support1800Tab$protein1,
                             support1800Tab$protein2),
                         static = edgeList(staticTab$protein1,
                             staticTab$protein2),
                         dpim1=edgeStrings$dpim1)

}

{
    vennFile = "/home/glocke/DPiM/plots/consensus5_vsStatic.sharedEdges.tiff"
    mn <- list(cons = "support1800", static="Static", dpim1="DPIM1")
    venn.diagram(edgeStrings2, filename=vennFile, sub.fontfamily="sans-serif",
                 main.fontfamily="sans-serif",
                 category.names = mn, fill=rainbow(3),
                 alpha=c(0.5, 0.5, 0.5))
}



{
    cons334Tab <- read.table("/home/glocke/DPiM/dpim4/withInstr/consensus5_nets/support0334/support0334.net", header = T)
    nrTab <- read.table("/home/glocke/DPiM/dpim4/withInstr/nrBait/nrBait.net", header = T)
    dpim1Tab <- read.table("/home/glocke/DPiM/prevDPIM/DPIM1_scores.updateFBgn.tsv", header = T)

    edgeStrings <- list(
        cons=edgeList(cons334Tab$protein1, cons334Tab$protein2),
        nr=edgeList(nrTab$protein1, nrTab$protein2),
        dpim1=edgeList(dpim1Tab$Interactor_1, dpim1Tab$Interactor_2)
        )
}

{
    vennFile = "/home/glocke/DPiM/plots/consensus5_nrBait.sharedEdges.tiff"
    mn <- list(cons="support334", nr="nrBait", dpim1="DPIM1")
    venn.diagram(edgeStrings, filename=vennFile, sub.fontfamily="sans-serif",
                 main.fontfamily="sans-serif",
                 category.names = mn, fill=rainbow(3),
                 alpha=c(0.5, 0.5, 0.5))
}

{
    nrTab <- read.table("/home/glocke/DPiM/dpim4/withInstr/nrBait/nrBait.net", header = T)
    dpim1Tab <- read.table("/home/glocke/DPiM/prevDPIM/DPIM1_scores.updateFBgn.tsv", header = T)
    newTab <- read.table("/home/glocke/DPiM/newDpim2/nrBait/newDpim2.nrBait.net", header = T)
    edgeStrings <- list(
        dpim1=edgeList(dpim1Tab$Interactor_1, dpim1Tab$Interactor_2),
        dpim1.2=edgeList(newTab$protein1, newTab$protein2),
        dpim4=edgeList(nrTab$protein1, nrTab$protein2)
        )
    for (n in names(edgeStrings)) {
        outFile <- paste("~/tmp/", n, ".edge.txt", sep="");
        write(edgeStrings[[n]], file=outFile)
    }
}

{
    vennFile = "/home/glocke/DPiM/plots/remappedDpim1.sharedEdges.tiff"
    mn <- list(dpim1="DPIM1", dpim1.2="DPIM1.2", dpim4="bestRep")
    venn.diagram(edgeStrings, filename=vennFile, sub.fontfamily="sans-serif",
                 main.fontfamily="sans-serif",
                 category.names = mn, fill=rainbow(3),
                 alpha=c(0.5, 0.5, 0.5))
}

{
    edge2 <- edgeStrings
    edge2$dpim4 <- NULL
    edge2 <- rev(edge2)
    vennFile = "/home/glocke/DPiM/plots/remappedDpim1only.sharedEdges.tiff"
    mn <- list(dpim1.2="DPIM1.2", dpim1="DPIM1")
    venn.diagram(edge2, filename=vennFile, sub.fontfamily="sans-serif",
                 main.fontfamily="sans-serif",
                 category.names = mn, fill=rainbow(2),
                 alpha=c(0.5, 0.5))
}

{
    dpim12Tab <- read.table("/home/glocke/DPiM/newDpim2/nrBait/newDpim2.nrBait.net", header = T)
    nrTab <- read.table("/home/glocke/DPiM/dpim4/withInstr/nrBait/nrBait.net", header = T)
    meanTab <- read.table("/home/glocke/DPiM/dpim4/withInstr/meanBait/meanBait.net", header = T)
    edgeStrings <- list(
        dpim1.2=edgeList(dpim12Tab$protein1, dpim12Tab$protein2),
        nr=edgeList(nrTab$protein1, nrTab$protein2),
        mean=edgeList(meanTab$protein1, meanTab$protein2)
        )
    
}

{
    vennFile = "/home/glocke/DPiM/plots/nrBait_vs_meanBait_1.2.sharedEdges.tiff"
    mn <- list(dpim1.2="DPIM1.2", nr="bestRep", mean="allRep")
    venn.diagram(edgeStrings, filename=vennFile, sub.fontfamily="sans-serif",
                 main.fontfamily="sans-serif",
                 category.names = mn, fill=rainbow(3),
                 alpha=c(0.5, 0.5, 0.5))
}

{
    vennFile = "/home/glocke/DPiM/plots/nrBait_meanBait_newDpim.sharedEdges.tiff"
    mn <- list(nr="bestRep", mean="allRep", dpim1="DPIM1", dpim1.2="DPIM1.2")
    edgeStrings <- list(edgeStrings$nr, edgeStrings$mean, edgeStrings$dpim1,
                        edgeStrings$dpim1.2)
    venn.diagram(edgeStrings, filename=vennFile, sub.fontfamily="sans-serif",
                 main.fontfamily="sans-serif",
                 category.names = mn, fill=rainbow(4),
                 alpha=0.5)
}


{
    hgScore <- read.table("/home/glocke/DPiM/human/biopStyleData/nrBait/biopNrBait.net", header = T)
    comp85 <- read.table("/home/glocke/DPiM/human/CompPASS/BioPlex_8467map-edges.net", header = T)
    comp59 <- read.table("/home/glocke/DPiM/human/CompPASS/BioPlex.v5884_HCIs.net", header = T)
    edgeStrings <- list(
        c59=edgeList(comp59$protein1, comp59$protein2),
        c85=edgeList(comp85$protein1, comp85$protein2),
        hg=edgeList(hgScore$protein1, hgScore$protein2)
        )
    
}

{
    vennFile = "/home/glocke/DPiM/plots/hgVsCompPASS.sharedEdges_08-23-2016.tiff"
    mn <- list(c59="BioPlex\n(5884)", c85="CompPASS\n(8467)", hg="Human-HG (8467)")
    venn.diagram(edgeStrings, filename=vennFile, sub.fontfamily="sans-serif",
                 main.fontfamily="sans-serif",
                 category.names = mn, fill=rainbow(3),
                 alpha=c(0.5, 0.5, 0.5))
}

{
    nodeStrings <- list(
        c59=unique(c(comp59$protein1, comp59$protein2)),
        c85=unique(c(comp85$protein1, comp85$protein2)),
        hg=unique(c(hgScore$protein1, hgScore$protein2))
        )
    vennFile = "/home/glocke/DPiM/plots/hgVsCompPASS.sharedNodes_08-23-2016.tiff"
    mn <- list(c59="BioPlex\n(5884)", c85="CompPASS\n(8467)", hg="Human-HG (8467)")
    venn.diagram(nodeStrings, filename=vennFile, sub.fontfamily="sans-serif",
                 main.fontfamily="sans-serif",
                 category.names = mn, fill=rainbow(3),
                 alpha=c(0.5, 0.5, 0.5))
    
}

{
    hgBB <- read.table("/home/glocke/DPiM/human/biopStyleData/nrBait/biopNrBait.baitBait", header = T)
    hgBP <- read.table("/home/glocke/DPiM/human/biopStyleData/nrBait/biopNrBait.baitPrey", header = T)
    hgPP <- read.table("/home/glocke/DPiM/human/biopStyleData/nrBait/biopNrBait.preyPrey", header = T)
    
    edgeStrings2 <- list(
        bb = edgeList(hgBB$protein1, hgBB$protein2),
        bp = edgeList(hgBP$protein1, hgBP$protein2),
        ##pp = edgeList(hgPP$protein1, hgPP$protein2)
        )
   edgeStrings2$c85 <- edgeStrings$c85
}

{
    vennFile = "/home/glocke/DPiM/plots/hgVsCompPASS.sharedEdgesByType_08-23-2016.tiff"
    mn <- list(bb="Human-HG bait-bait", bp="Human-HG bait-prey",
               c85="CompPASS\n(8467)")
    venn.diagram(edgeStrings2, filename=vennFile, sub.fontfamily="sans-serif",
                 main.fontfamily="sans-serif",
                 category.names = mn, fill=rainbow(3),
                 alpha=c(0.5, 0.5, 0.5))
}

{
    dpim1Tab <- read.table("/home/glocke/DPiM/prevDPIM/DPIM1_scores.updateFBgn.tsv", header = T)
    dpim12Tab <- read.table("/home/glocke/DPiM/newDpim1_2/nrBait/newDpim2.nrBait.net", header = T)
    dpim2Tab <- read.table("/home/glocke/DPiM/prevDPIM/dpim2_nrtap.120123.nrBait.58.17.network", header = T)
    dpim22Tab <- read.table("/home/glocke/DPiM/newDpim2/nrBait_08-23-2016/newDpim2_08-2016.net", header = T)
    edgeStrings <- list(
        dpim1=edgeList(dpim1Tab$Interactor_1, dpim1Tab$Interactor_2),
        dpim1.2=edgeList(dpim12Tab$protein1, dpim12Tab$protein2),
        dpim2=edgeList(dpim2Tab$node1, dpim2Tab$node2),
        dpim2.2=edgeList(dpim22Tab$protein1, dpim22Tab$protein2),
        )
    
}

{
    vennFile = "/home/glocke/DPiM/plots/remappedDPiM1and2.sharedEdges.tiff"
    mn <- list(dpim1="DPIM1", dpim1.2="DPIM1.2", dpim2="DPIM2", dpim2.2="DPIM2.2")
    venn.diagram(edgeStrings, filename=vennFile, sub.fontfamily="sans-serif",
                 main.fontfamily="sans-serif",
                 category.names = mn, fill=rainbow(4),
                 alpha=0.5)
}


{
    dpim1Tab <- read.table("/home/glocke/DPiM/prevDPIM/DPIM1_scores.updateFBgn.tsv", header = T)
    dpim1JulTab <- read.table("/home/glocke/DPiM/newDpim1_2/newNrBait/newNrBait.net", header = T)
    dpim1Aug1Tab <- read.table("/home/glocke/DPiM/augRemap/apmsData/dpim1/nrBait09-12-2016/nrBait.net", header = T)
    dpim1Aug0Tab <- read.table("/home/glocke/DPiM/augRemap/apmsData/dpim1/oldNrBait09-12-2016/oldNrBait.net", header = T)
    edgeStrings <- list(
        dpim1=edgeList(dpim1Tab$Interactor_1, dpim1Tab$Interactor_2),
        dpim1Jul=edgeList(dpim1JulTab$protein1, dpim1JulTab$protein2),
        dpim1Aug1=edgeList(dpim1Aug1Tab$protein1, dpim1Aug1Tab$protein2),
        dpim1Aug0=edgeList(dpim1Aug0Tab$protein1, dpim1Aug0Tab$protein2)
        )
    
}

{
    vennFile = "/home/glocke/DPiM/plots/remappedDPiM1versions.sharedEdges.tiff"
    mn <- list(dpim1="Published", dpim1Jul="July 2016", dpim1Aug1="Current\npipeline", dpim1Aug0="Previous\npipeline")
    venn.diagram(edgeStrings, filename=vennFile, sub.fontfamily="sans-serif",
                 main.fontfamily="sans-serif",
                 category.names = mn, fill=rainbow(4),
                 alpha=0.5)
}


{
    dpim1Tab <- read.table("/home/glocke/DPiM/prevDPIM/DPIM1_scores.updateFBgn.tsv", header = T)
    dpim1NincTab <- read.table("/home/glocke/DPiM/augRemap/apmsData/dpim1/newNrBait09-13-2016/newNrBait.net", header = T)
    dpim1Aug1Tab <- read.table("/home/glocke/DPiM/augRemap/apmsData/dpim1/nrBait09-12-2016/nrBait.net", header = T)
    edgeStrings <- list(
        dpim1=edgeList(dpim1Tab$Interactor_1, dpim1Tab$Interactor_2),
        dpim1Ninc=edgeList(dpim1NincTab$protein1, dpim1NincTab$protein2),
        dpim1Dec=edgeList(dpim1Aug1Tab$protein1, dpim1Aug1Tab$protein2)
        )
    
}

{
    vennFile = "/home/glocke/DPiM/plots/remappedDPiM1_compareLC.sharedEdges.tiff"
    mn <- list(dpim1="Published", dpim1Ninc="Non-increasing", dpim1Aug0="Decreasing")
    venn.diagram(edgeStrings, filename=vennFile, sub.fontfamily="sans-serif",
                 main.fontfamily="sans-serif",
                 category.names = mn, fill=rainbow(3),
                 alpha=0.5)
}
