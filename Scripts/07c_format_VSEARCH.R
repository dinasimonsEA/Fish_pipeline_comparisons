# script to format the VSEARCH outputs so we can assign taxonomy

# get path
path <- getwd()

# Dataset-specific directories
output_dir   <- file.path(path, "Data", "Temp", test_data_name)
rds_dir      <- file.path(output_dir, "R_objects")
vsearch_dir  <- file.path(output_dir, "06b_vsearch_denovo_outputs", "05_nochim")
processed_dir <- file.path(path, "Data", "Processed", test_data_name)

# ensure dirs exist
dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)

# get VSEARCH fasta files
files <- list.files(vsearch_dir, pattern = "\\.fasta$", full.names = TRUE)

if (length(files) == 0) {
  stop("No VSEARCH nochim fasta files found — check previous steps.")
}

# function to read one fasta and return named vector (sequence -> count)
read_vsearch_fasta <- function(file) {
  dna <- readDNAStringSet(file)
  headers <- names(dna)

  # extract counts from "size=XXXX"
  counts <- as.numeric(sub(".*size=([0-9]+).*", "\\1", headers))
  seqs <- as.character(dna)

  # combine identical sequences just in case
  tapply(counts, seqs, sum)
}

# read all samples
sample_list <- lapply(files, read_vsearch_fasta)

names(sample_list) <- sub(
  "_merged_filtered_uniques_denoised_nochim\\.fasta$",
  "",
  basename(files)
)

# get all unique sequences (ASVs)
all_seqs <- unique(unlist(lapply(sample_list, names)))

# build count matrix (samples x ASVs)
seqtab <- matrix(
  0,
  nrow = length(sample_list),
  ncol = length(all_seqs),
  dimnames = list(names(sample_list), all_seqs)
)

# fill matrix
for (i in seq_along(sample_list)) {
  seqs <- names(sample_list[[i]])
  counts <- sample_list[[i]]
  seqtab[i, seqs] <- counts
}

# save R object
saveRDS(seqtab, file = file.path(rds_dir, "06_seqtab_VSEARCH.rds"))

# ==============================
# MATCH DADA2 OUTPUT FORMAT
# ==============================

# FASTA output
sequences <- colnames(seqtab)
headers <- paste0(">ASV_", seq_along(sequences))
output_fasta <- c(rbind(headers, sequences))

write(
  output_fasta,
  file = file.path(processed_dir, "06b_ASV_seqs_VSEARCH.fasta")
)

# TSV counts
tab <- t(seqtab)

write.table(
  tab,
  file = file.path(processed_dir, "06b_ASV_counts_VSEARCH.tsv"),
  sep = "\t",
  quote = FALSE,
  col.names = NA
)