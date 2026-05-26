# script for dereplicating, merging and chimera removal using DADA2

# get path
path <- getwd()

# Dataset-specific directories
output_dir  <- file.path(path, "Data", "Temp", test_data_name)
rds_dir     <- file.path(output_dir, "R_objects")
processed_dir <- file.path(path, "Data", "Processed", test_data_name)

# ensure directories exist
dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)

# read in filtered files
filtFs <- readRDS(file = file.path(rds_dir, "04_fnFs.filtN.rds"))
filtRs <- readRDS(file = file.path(rds_dir, "04_fnRs.filtN.rds"))

# read in error models
errF <- readRDS(file = file.path(rds_dir, "05_errF.rds"))
errR <- readRDS(file = file.path(rds_dir, "05_errR.rds"))

# check files exist
exists <- file.exists(filtFs) & file.exists(filtRs)

if (!any(exists)) {
  stop("No valid filtered files available for DADA2 step.")
}

# dereplicate the reads
derepFs <- derepFastq(filtFs[exists])
derepRs <- derepFastq(filtRs[exists])

# get sample names
sample.names <- readRDS(file = file.path(rds_dir, "03_sample_names.rds"))

# name the derep objects
names(derepFs) <- sample.names[exists]
names(derepRs) <- sample.names[exists]

# perform DADA2 sample inference
dadaFs <- dada(derepFs, err = errF, multithread = TRUE)
dadaRs <- dada(derepRs, err = errR, multithread = TRUE)

# merge paired reads
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose = TRUE)

# make ASV table
seqtab <- makeSequenceTable(mergers)

# remove chimeras
seqtab.nochim <- removeBimeraDenovo(
  seqtab,
  method = "consensus",
  multithread = TRUE,
  verbose = TRUE
)

# write ASV sequences (FASTA)
sequences <- colnames(seqtab.nochim)
headers <- paste0(">ASV_", seq_along(sequences))

output_fasta <- c(rbind(headers, sequences))

write(
  output_fasta,
  file = file.path(processed_dir, "06_ASV_seqs_DADA2.fasta")
)

# write ASV counts
tab <- t(seqtab.nochim)

write.table(
  tab,
  file = file.path(processed_dir, "06_ASV_counts_DADA2.tsv"),
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

# save R objects
saveRDS(dadaFs, file = file.path(rds_dir, "06_dadaFs.rds"))
saveRDS(dadaRs, file = file.path(rds_dir, "06_dadaRs.rds"))
saveRDS(mergers, file = file.path(rds_dir, "06_mergers.rds"))
saveRDS(seqtab.nochim, file = file.path(rds_dir, "06_seqtab.nochim.rds"))