## if there are errors when reading the file, the input file needs to be
## manually corrected to make sure each line starts with "w"
argv <- commandArgs(TRUE)
inFile <- argv[1]
outFile <- argv[2]
#inFile <- '/home/glocke/DPiM/dpim3_raw/dpim2_all.120123.plusUniqMsInst.cp'
test <- read.table(inFile, sep="\t", header=T,
                   quote="", comment.char="",stringsAsFactors=F)

## date string would be missing a 0 if the month is 1-9.
test$date1 <- paste0("0", test$date[nchar(test$date)==5])
test$date <- as.Date(test$date1,"%m%d%y")

## replace empty bait rows
test$bait[test$bait==""] = 'x_FBgn0000000'

## replace bait rows with 2 "_"'s
test$bait[grep("_.*_",test$bait)] = 'x_GFP'

## new data frame
newDF <- data.frame("tap_id" = test$tap_id,
                    "ms_inst_run_id" = test$name,
                    "user"="Automatic",
                    "search_id" = test$search_id,
                    sample_date = test$date,
                    "total_peptides" = test$abundance,
                    "unique_peptides" = "NaN",
                    "bait_ref" = matrix(unlist(strsplit(test$bait,"_")),
                        ncol=2,byrow=T)[,2],
                    "prey_ref" = test$gene_symbol)
                    

## write out the data frame
write.table(newDF, file=outFile, sep="\t", row.names=F, quote=F)
            

