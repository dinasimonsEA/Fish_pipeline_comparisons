# assigning taxonomy using dada2 and metabeat rdp
# vsearch in other script

# get path
path <- getwd()

# Dataset-specific directories
output_dir   <- file.path(path, "Data", "Temp", test_data_name)
rds_dir      <- file.path(output_dir, "R_objects")
processed_dir <- file.path(path, "Data", "Processed", test_data_name)

# ensure directories exist
dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)

# read in the R object
seqtab.nochim <- readRDS(file = file.path(rds_dir, "06_seqtab.nochim.rds"))

# database path
database <- file.path(path, "Data", "Databases", "EA_custom", "EA_fish_database.riaz.AssignTaxonomy.fasta")

# jono to fix nBoot, should be 1000
# assign taxonomy (dada2)
taxa <- assignTaxonomy(
  seqtab.nochim,
  database,
  tryRC = TRUE,
  verbose = TRUE,
  multithread = TRUE,
  minBoot = 90,
  outputBootstraps = TRUE
)

boot_dada2<-taxa[["boot"]]
taxa_dada2<-taxa[["tax"]]
taxa_dada2_print <- taxa_dada2
rownames(taxa_dada2_print) <- NULL

# database name (used in filenames)
db_name <- "EA_fish_riaz"

# write taxonomy table
write.table(
  taxa_dada2,
  file = file.path(processed_dir, paste0("08_assigned_taxonomy_DADA2_", db_name, ".csv")),
  sep = ",",
  quote = FALSE
)

write.table(
  boot_dada2,
  file = file.path(processed_dir, paste0("08_bootstraps_DADA2_", db_name, ".csv")),
  sep = ",",
  quote = FALSE
)

# save R object
saveRDS(taxa_dada2, file = file.path(rds_dir, "08_taxa.rds"))