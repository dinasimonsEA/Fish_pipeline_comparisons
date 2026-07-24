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
BLAST_MetaBEAT_output <- read.table(
  file = file.path(processed_dir, "08_ASVs_MetaBEAT_blast_EA_riaz.txt"),
  sep = "\t",
  header = FALSE,
  stringsAsFactors = FALSE
)

# clean species column
BLAST_MetaBEAT_output <- BLAST_MetaBEAT_output %>%
  mutate(
    V3 = sub("^[^ ]+ ", "", V3),  # remove accession prefix
    V3 = gsub("_", " ", V3)       # optional: make names readable
  )

# set column names (BLAST outfmt 6)
colnames(BLAST_MetaBEAT_output) <- c(
  "qseqid", "sseqid", "species", "pident", "length",
  "mismatch", "gapopen", "qstart", "qend",
  "sstart", "send", "evalue", "bitscore"
)

# join taxonomy lookup
BLAST_MetaBEAT_output <- BLAST_MetaBEAT_output %>%
  left_join(taxonomy_lookup, by = c("sseqid" = "gbAccession"))

# 4. Filter ---------
BLAST_MetaBEAT_output <- BLAST_MetaBEAT_output %>%
  filter(pident >= min_pident) %>%
  filter(length >= min_length) %>%
  group_by(qseqid) %>%
  mutate(max_bitscore = max(bitscore)) %>%
  filter(bitscore >= (top_frac * max_bitscore)) %>%
  ungroup()

# 5. Get taxonomy matrix
## MetaBEAT
MetaBEAT_taxa <- BLAST_MetaBEAT_output %>%
select(c(kingdom, phylum, class, order, family, genus, species))
colnames(MetaBEAT_taxa) <- 
 c(
  "Kingdom",
  "Phylum",
  "Class",
  "Order",
  "Family",
  "Genus",
  "Species"
)

MetaBEAT_taxa_matrix <- as.matrix(MetaBEAT_taxa)
head(MetaBEAT_taxa_matrix)

# 6. Run LCA
combined_data_MetaBEAT <- cbind(ASV_ID = BLAST_MetaBEAT_output$qseqid, as.data.frame(MetaBEAT_taxa))

lca_results_MetaBEAT <- combined_data_MetaBEAT %>%
  group_by(ASV_ID) %>%
  do({
    #Use distinct() here to save massive amounts of time/temp space
    sub_matrix <- .[, -1] %>% distinct() %>% as.matrix()
    lca <- taxonomizr::condenseTaxa(sub_matrix)
    as.data.frame(lca)
  })

# 7. Final formatting
## MetaBEAT
lca_results_MetaBEAT[lca_results_MetaBEAT == "NA"] <- NA

lca_results_MetaBEAT <- lca_results_MetaBEAT %>%
  ungroup() %>%
  mutate(
    lca_rank = apply(.[, -1], 1, function(x) {
      idx <- which(!is.na(x) & x != "")
      if (length(idx) > 0) colnames(lca_results_MetaBEAT)[-1][max(idx)] else NA
    }),
    lca_name = apply(.[, -1], 1, function(x) {
      idx <- which(!is.na(x) & x != "")
      if (length(idx) > 0) x[max(idx)] else NA
    })
  ) %>%
  select(ASV_ID, everything(), lca_rank, lca_name)

# 7. Save Outputs
out_file_MetaBEAT <- paste0(processed_dir, "/09_MetaBEAT_ASV_LCA_results_EA_riaz", min_pident, "_", min_length, "_top_", top_frac, "%.tsv")
write.table(lca_results_MetaBEAT, file = out_file_MetaBEAT, sep = "\t", quote = FALSE, row.names = FALSE)
