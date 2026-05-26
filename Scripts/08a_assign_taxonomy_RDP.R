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

# database path
database <- file.path(
  path,
  "Data", "Databases", "Meta-fish-lib",
  "references.12s.miya.dada.taxonomy.v268.fasta"
)

# assign taxonomy
taxa <- assignTaxonomy(
  seqtab.nochim,
  database,
  tryRC = TRUE,
  verbose = TRUE,
  multithread = TRUE
)

# database name (used in filenames)
db_name <- "Meta-fish-lib"

# write taxonomy table
write.table(
  taxa,
  file = file.path(processed_dir, paste0("08_assigned_taxonomy_DADA2_", db_name, ".csv")),
  sep = ",",
  quote = FALSE
)

# save R object
saveRDS(taxa, file = file.path(rds_dir, "08_taxa.rds"))