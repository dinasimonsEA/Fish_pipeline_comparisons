# script to format the VSEARCH outputs so we can assign taxonomy

# folder containing all VSEARCH fasta files (one per sample)
input_dir <- "Data/Temp/06b_vsearch_denovo_outputs/05_nochim"
files <- list.files(input_dir, pattern="\\.fasta$", full.names=TRUE)

# function to read one fasta and return named vector (sequence -> count)
read_vsearch_fasta <- function(file) {
  dna <- readDNAStringSet(file)
  headers <- names(dna)

  # extract counts from "size=XXXX"
  counts <- as.numeric(sub(".*size=([0-9]+)", "\\1", headers))
  seqs <- as.character(dna)

  # combine identical sequences just in case
  tapply(counts, seqs, sum)
}

# read all samples
sample_list <- lapply(files, read_vsearch_fasta)
names(sample_list) <- sub("_merged_filtered_uniques_denoised_nochim\\.fasta$", "", basename(files))

# get all unique sequences (ASVs)
all_seqs <- unique(unlist(lapply(sample_list, names)))

# build count matrix (samples x ASVs)
seqtab <- matrix(0,
                 nrow = length(sample_list),
                 ncol = length(all_seqs),
                 dimnames = list(names(sample_list), all_seqs))

# fill matrix
for (i in seq_along(sample_list)) {
  seqs <- names(sample_list[[i]])
  counts <- sample_list[[i]]
  seqtab[i, seqs] <- counts
}

saveRDS(seqtab, file = paste(path, "/Data/Temp/R_objects/06_seqtab_VSEARCH.rds", sep=""))

# MATCH DADA2 OUTPUT

# FASTA output
sequences <- colnames(seqtab)
headers <- paste0(">ASV_", seq_along(sequences))
output_fasta <- c(rbind(headers, sequences))

write(output_fasta, file = "Data/Processed/06b_ASV_seqs_VSEARCH.fasta")

# TSV counts
tab <- t(seqtab)
write.table(tab, file = "Data/Processed/06b_ASV_counts_VSEARCH.tsv", sep = "\t", quote = FALSE, col.names = NA)