# Script to filter the data to remove any poor quality reads

# get wd & data location
path<-getwd()
test_data_loc <- paste("Data/Raw/RingTrial_Sean")

# read in the lists of cutadapt files for trimming
fnFs.cut <- readRDS(file = paste(path,"/Data/Temp/R_objects/02_fnFs.cut.rds",sep=""))
fnRs.cut <- readRDS(file = paste(path,"/Data/Temp/R_objects/02_fnRs.cut.rds",sep=""))

# set up filepaths for output files
fnFs.filtN <- file.path(path, "/Data/Temp/filterAndTrim", basename(fnFs.cut))
fnRs.filtN <- file.path(path, "/Data/Temp/filterAndTrim", basename(fnRs.cut))

# run filterAndTrim based on user defined parameters
out <- filterAndTrim(fnFs.cut, fnFs.filtN, fnRs.cut, fnRs.filtN,
maxN = 0, maxEE = maxEE, truncQ = truncQ, minLen = minLen, truncLen = truncLen, multithread = TRUE) 

# write objects to pass to next script
saveRDS(fnFs.filtN, file = paste(path, "/Data/Temp/R_objects/04_fnFs.filtN.rds" ,sep=""))
saveRDS(fnRs.filtN, file = paste(path, "/Data/Temp/R_objects/04_fnRs.filtN.rds" ,sep=""))
saveRDS(out, file = paste(path, "/Data/Temp/R_objects/04_out.rds", sep=""))

# generate quality plots on the quality-trimmed data
pdf(file = paste(path,"/Results/04_post_trim_quality_plots.pdf", sep=""), 
        width = 10,
        height = 10)