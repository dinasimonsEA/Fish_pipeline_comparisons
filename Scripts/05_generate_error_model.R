# Script to generate error models

# get wd & data location
path<-getwd()
test_data_loc <- paste("Data/Raw/RingTrial_Sean")

# read in the lists of filtered reads
filtFs <- readRDS(file = paste(path, "/Data/Temp/R_objects/04_fnFs.filtN.rds" ,sep=""))
filtRs <- readRDS(file = paste(path, "/Data/Temp/R_objects/04_fnRs.filtN.rds" ,sep=""))
out <- readRDS(file = paste(path, "/Data/Temp/R_objects/04_out.rds", sep=""))

# check that the files in filtFs and filtRs lists exist - I.e. made it through filterAndTrim
# (sometimes files can have no reads pass and it upsets the next stage)
exists <- file.exists(filtFs)

# run error rates but only on those files that exist

## learn the error rates
errF <- learnErrors(filtFs[exists], multithread = TRUE, nbases = 2e+08,)
errR <- learnErrors(filtRs[exists], multithread = TRUE, nbases = 2e+08)

## write out error rates for use later
saveRDS(errF, file = paste(path,"/Data/Temp/R_objects/05_errF.rds",sep=""))
saveRDS(errR, file = paste(path,"/Data/Temp/R_objects/05_errR.rds",sep=""))

# write plot to file for inspection
pdf(file = paste(path,"/Results/05_error_rate_plots.pdf", sep=""), width = 10, height = 10)