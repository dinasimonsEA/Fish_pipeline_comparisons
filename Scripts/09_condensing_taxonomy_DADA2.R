# script for condensing taxonomy for BLAST outputs using LCA
# uses taxonomiser

# 1. File paths ---------
path <- getwd()
processed_dir <- file.path(path, "Data", "Processed", test_data_name)

# 2. Make taxonomy & accession look-up table from metabeat ---------
## read MetaFishLib database
mfl_db <- read.csv(
  file = file.path(path, "Data", "Databases", "EA_custom", "EA_fish_database.riaz.csv"),
  stringsAsFactors = FALSE
)

## make taxonomy lookup table
taxonomy_lookup <- mfl_db %>%
  select(gbAccession, genus, family, order, class, phylum, kingdom)

# 3. Load blast data and tidy ---------
## DADA2 BLAST output
BLAST_DADA2_output <- read.table(
  file = file.path(processed_dir, "08_ASVs_DADA2_blast_EA_riaz.txt"),
  sep = "\t",
  header = FALSE,
  stringsAsFactors = FALSE
)

# clean species column
BLAST_DADA2_output <- BLAST_DADA2_output %>%
  mutate(
    V3 = sub("^[^ ]+ ", "", V3),  # remove accession prefix
    V3 = gsub("_", " ", V3)       # optional: make names readable
  )

# set column names (BLAST outfmt 6)
colnames(BLAST_DADA2_output) <- c(
  "qseqid", "sseqid", "species", "pident", "length",
  "mismatch", "gapopen", "qstart", "qend",
  "sstart", "send", "evalue", "bitscore"
)

# join taxonomy lookup
BLAST_DADA2_output <- BLAST_DADA2_output %>%
  left_join(taxonomy_lookup, by = c("sseqid" = "gbAccession"))

# 4. Filter ---------
BLAST_DADA2_output <- BLAST_DADA2_output %>%
  filter(pident >= min_pident) %>%
  filter(length >= min_length) %>%
  group_by(qseqid) %>%
  mutate(max_bitscore = max(bitscore)) %>%
  filter(bitscore >= (top_frac * max_bitscore)) %>%
  ungroup()

# 5. Get taxonomy matrix
dada2_taxa <- BLAST_DADA2_output %>%
select(c(kingdom, phylum, class, order, family, genus, species))
colnames(dada2_taxa) <- 
 c(
  "Kingdom",
  "Phylum",
  "Class",
  "Order",
  "Family",
  "Genus",
  "Species"
)

dada2_taxa_matrix <- as.matrix(dada2_taxa)
head(dada2_taxa_matrix)

# 6. Run LCA
combined_data_dada2 <- cbind(ASV_ID = BLAST_DADA2_output$qseqid, as.data.frame(dada2_taxa))

lca_results_dada2 <- combined_data_dada2 %>%
  group_by(ASV_ID) %>%
  do({
    # Use distinct() here to save massive amounts of time/temp space
    sub_matrix <- .[, -1] %>% distinct() %>% as.matrix()
    lca <- taxonomizr::condenseTaxa(sub_matrix)
    as.data.frame(lca)
  })

# 7. Final formatting
lca_results_dada2[lca_results_dada2 == "NA"] <- NA

lca_results_dada2 <- lca_results_dada2 %>%
  ungroup() %>%
  mutate(
    lca_rank = apply(.[, -1], 1, function(x) {
      idx <- which(!is.na(x) & x != "")
      if (length(idx) > 0) colnames(lca_results_dada2)[-1][max(idx)] else NA
    }),
    lca_name = apply(.[, -1], 1, function(x) {
      idx <- which(!is.na(x) & x != "")
      if (length(idx) > 0) x[max(idx)] else NA
    })
  ) %>%
  select(ASV_ID, everything(), lca_rank, lca_name)

# 7. Save Outputs
out_file_dada2 <- paste0(processed_dir, "/09_DADA2_ASV_LCA_results_EA_riaz", min_pident, "_", min_length, "_top_", top_frac, "%.tsv")
write.table(lca_results_dada2, file = out_file_dada2, sep = "\t", quote = FALSE, row.names = FALSE)