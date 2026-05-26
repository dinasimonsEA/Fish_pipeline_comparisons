# script for sequence tracking dada2

# get path
path <- getwd()

# Dataset-specific directories
output_dir  <- file.path(path, "Data", "Temp", test_data_name)
rds_dir     <- file.path(output_dir, "R_objects")
results_dir <- file.path(path, "Results", test_data_name)

# ensure results dir exists
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

# read in the R objects
out            <- readRDS(file = file.path(rds_dir, "04_out.rds"))
dadaFs         <- readRDS(file = file.path(rds_dir, "06_dadaFs.rds"))
dadaRs         <- readRDS(file = file.path(rds_dir, "06_dadaRs.rds"))
mergers        <- readRDS(file = file.path(rds_dir, "06_mergers.rds"))
seqtab.nochim  <- readRDS(file = file.path(rds_dir, "06_seqtab.nochim.rds"))
filtFs         <- readRDS(file = file.path(rds_dir, "04_fnFs.filtN.rds"))

# ensure only existing files are used (important!)
exists <- file.exists(filtFs)

# build tracking table
track <- cbind(
  out[exists, 1],
  out[exists, 2],
  sapply(dadaFs, getN),
  sapply(dadaRs, getN),
  sapply(mergers, getN),
  rowSums(seqtab.nochim)
)

colnames(track) <- c(
  "input", 
  "filtered",
  "denoisedF",
  "denoisedR",
  "merged",
  "nochim"
)

# get sample names (aligned to filtered data)
sample.names <- manifest$sampleID[exists]
rownames(track) <- sample.names

# write tracking table
write.csv(
  track,
  file = file.path(results_dir, "07_track_reads_table.csv"),
  row.names = TRUE
)

print(track)