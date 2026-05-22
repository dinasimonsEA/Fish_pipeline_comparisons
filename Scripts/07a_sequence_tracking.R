# script for sequence tracking dada2

# get path
path<-getwd()

# read in the R objects
out <- readRDS(file = paste(path, "/Data/Temp/R_objects/04_out.rds", sep=""))
dadaFs <- readRDS(file = paste(path, "/Data/Temp/R_objects/06_dadaFs.rds", sep=""))
dadaRs <- readRDS(file = paste(path, "/Data/Temp/R_objects/06_dadaRs.rds", sep=""))
mergers	<- readRDS(file = paste(path, "/Data/Temp/R_objects/06_mergers.rds", sep=""))
seqtab.nochim <- readRDS(file = paste(path, "/Data/Temp/R_objects/06_seqtab.nochim.rds", sep=""))
filtFs <- readRDS(file = paste(path, "/Data/Temp/R_objects/04_fnFs.filtN.rds", sep=""))

track <- cbind(
	out[,1],
	out[,2],
	sapply(dadaFs, getN),
	sapply(dadaRs, getN),
	sapply(mergers, getN),
	rowSums(seqtab.nochim))
colnames(track) <- c(
	"input", 
	"filtered",
	"denoisedF",
	"denoisedR",
	"merged",
	"nochim")

# get sample names
sample.names <- manifest$sampleID
rownames(track) <- sample.names

# write tracking table
write.table(track, file = paste(path, "/Results/07_track_reads_table.csv", sep=""))
print(track)