# get path
path<-getwd()

# read in the R objects
seqtab.nochim <- readRDS(file = paste(path, "/Data/Temp/R_objects/06_seqtab.nochim.rds", sep=""))

# find path for db
database <- "Data/Databases/Meta-fish-lib/references.12s.miya.dada.taxonomy.v268.fasta"

# assign taxonomy
taxa <- assignTaxonomy(seqtab.nochim, database, tryRC = TRUE, verbose = TRUE, multithread = TRUE)

## get db_name to add onto output
db_name <- "Meta-fish-lib"

# write taxonomy table
write.table(taxa, file = paste(path, "/Data/Processed/08_assigned_taxonomy_DADA2_", db_name, ".csv", sep=""))

# write out R objects for use later
saveRDS(taxa, file = paste(path, "/Data/Temp/R_objects/08_taxa.rds", sep=""))