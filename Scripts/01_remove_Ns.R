# Remove N's for DADA2 functionality using filter and trim

#get wd & data location
path<-getwd()
test_data_loc <- paste("Data/Raw/RingTrial_Sean")

# group raw files into "forward" reads and "reverse" reads
fnFs <- sort(manifest$absolute_forward_path)
fnRs <- sort(manifest$absolute_backward_path)

# set up filepaths for output files
fnFs.filtN <- file.path(path, "Data/Temp/01_filtN", basename(fnFs)) # Put N-filterd files in filtN/ subdirectory
fnRs.filtN <- file.path(path, "Data/Temp/01_filtN", basename(fnRs))

# run filterAndTrim for all files
eval(parse(text = paste("filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, ", ")")))

# ensure r object directory exists
dir.create(file.path(path, "Data/Temp/R_objects"), recursive = TRUE, showWarnings = FALSE)

# write objects to pass to next script
saveRDS(fnFs.filtN, file = paste(path,"/Data/Temp/R_objects/01_fnFs.filtN.rds",sep=""))
saveRDS(fnRs.filtN, file = paste(path,"/Data/Temp/R_objects/01_fnRs.filtN.rds",sep=""))