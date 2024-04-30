argv <- commandArgs(TRUE)
if (length(argv) != 1) {
    stop("usage: checkFDR.R <inFile>");
}
inFile <- argv[1]

test <- read.table(inFile, sep="\t", header=T, quote='',
                   stringsAsFactors=F)

print(paste("FDR",
            length(grep("reverse_",unique(c(test$prey_ref,test$bait_ref)))) /
            (length(unique(c(test$prey_ref,test$bait_ref))) -
                 length(grep("reverse_",unique(c(test$prey_ref,test$bait_ref))))),
            " "))
