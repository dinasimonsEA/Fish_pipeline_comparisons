# Script to explore quality of sequences 

# get wd & data location
path<-getwd()
test_data_loc <- paste("Data/Raw/RingTrial_Sean")

## read in the lists of cutadapt files
fnFs.cut <- readRDS(file = paste(path,"/Data/Temp/R_objects/02_fnFs.cut.rds",sep=""))
fnRs.cut <- readRDS(file = paste(path,"/Data/Temp/R_objects/02_fnRs.cut.rds",sep=""))

# extract sample names and write to R object
sample.names <- manifest$sampleID
print("Sample names:")
print(sample.names)
saveRDS(sample.names, file=paste(path, "/Data/Temp/R_objects/03_sample_names.rds", sep=""))

# write plot to file for inspection
pdf(file = paste(path,"/Results/03_pre_trim_quality_plots.pdf", sep=""),   # The directory you want to save the file into
    width = 10, # The width of the plot in inches
    height = 10)