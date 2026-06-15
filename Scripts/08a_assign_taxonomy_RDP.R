# assigning taxonomy using dada2 rdp

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

# read and format in METAbeat
seqtab.metabeat <- read.fasta(
  file = file.path(processed_dir, paste0("06_ASV_seqs_MetaBEAT.fasta"))
)

seqs_metabeat <- sapply(seqtab.metabeat, function(x) {
  paste(x, collapse = "")
})
seqs_metabeat <- toupper(seqs_metabeat)

# database path
database <- file.path(path, "Data", "Databases", "EA_custom", "EA_fish_database.riaz.AssignTaxonomy.fasta")

# assign taxonomy
taxa <- assignTaxonomy(
  seqtab.nochim,
  database,
  tryRC = TRUE,
  verbose = TRUE,
  multithread = TRUE
)

# assign taxonomy
taxa_metabeat <- assignTaxonomy(
  seqs_metabeat,
  database,
  tryRC = TRUE,
  verbose = TRUE,
  multithread = TRUE
)

# database name (used in filenames)
db_name <- "EA_fish_riaz"

# write taxonomy table
write.table(
  taxa,
  file = file.path(processed_dir, paste0("08_assigned_taxonomy_DADA2_", db_name, ".csv")),
  sep = ",",
  quote = FALSE
)

# write taxonomy table
write.table(
  taxa_metabeat,
  file = file.path(processed_dir, paste0("08_assigned_taxonomy_MetaBEAT_", db_name, ".csv")),
  sep = ",",
  quote = FALSE
)

# save R object
saveRDS(taxa, file = file.path(rds_dir, "08_taxa.rds"))
saveRDS(taxa_metabeat, file = file.path(rds_dir, "08_taxa_metabeat.rds"))