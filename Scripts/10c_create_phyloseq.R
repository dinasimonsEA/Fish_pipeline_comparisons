# script to take tidied long data and turn into phyloseq object

# ----------------------------------------------------------------
# create tax_table()

# =========================
# STEP 1: STRIP ALL Rle STRUCTURES
# =========================

taxa <- master_long_df_filtered %>%
  select(
    AphiaID_worms,
    taxa_name_final,
    kingdom,
    phylum,
    class,
    order,
    family,
    genus,
    species
  )

# FORCE everything to base types
taxa <- as.data.frame(lapply(taxa, function(x) {
  if (inherits(x, "Rle")) {
    return(as.vector(x))
  } else if (is.list(x)) {
    return(unlist(x))
  } else {
    return(x)
  }
}), stringsAsFactors = FALSE)

# convert back to df
taxa <- as.data.frame(taxa)

# =========================
# STEP 2: PRIORITY LOGIC
# =========================

taxa <- taxa %>%
  mutate(
    priority = case_when(
      !is.na(AphiaID_worms) & !is.na(kingdom) ~ 3,
      !is.na(kingdom) ~ 2,
      TRUE ~ 1
    )
  )

# =========================
# STEP 3: REMOVE DUPLICATES (SAFE)
# =========================

taxa <- taxa %>%
  group_by(taxa_name_final) %>%
  arrange(desc(priority), .by_group = TRUE) %>%
  slice_head(n = 1) %>%   # safer than slice(1)
  ungroup() %>%
  select(-priority)

# =========================
# STEP 4: FINAL CLEANUP
# =========================

taxa <- taxa %>%
  distinct() %>%
  as.data.frame()

rownames(taxa) <- NULL

# =========================
# STEP 5: SORT FOR PHYLOSEQ
# =========================

taxa$taxa_name_final <- as.factor(taxa$taxa_name_final)
taxa <- taxa[order(taxa$taxa_name_final), ]

# =========================
# STEP 6: TAX MATRIX
# =========================

tax_matrix <- as.data.frame(taxa)

rownames(tax_matrix) <- tax_matrix$taxa_name_final
tax_matrix$taxa_name_final <- NULL

tax_matrix <- as.matrix(tax_matrix)

# =========================
# CHECKS
# =========================

cat("Duplicate taxa:", any(duplicated(rownames(tax_matrix))), "\n")
cat("Missing kingdom:", sum(is.na(taxa$kingdom)), "\n")

# turn into phyloseq format
taxa <- as.data.frame(taxa)
taxa$taxa_name_final <- as.character(taxa$taxa_name_final)
taxa <- taxa[order(taxa$taxa_name_final), ]
rownames(taxa) <- taxa$taxa_name_final
taxa$taxa_name_final <- NULL
taxa <- as.matrix(taxa)

#-------------------------------------------------------------------------------
# create sample_data()
meta <- unique(subset(master_long_df_filtered, select = c("sampleid",
                                             "denoise_method",
                                             "taxonomy_method",
                                             "dataset",
                                             "fullID"
)))

meta <- as.data.frame(meta)
meta <- meta[order(meta$fullID), ]
rownames(meta) <- meta$fullID
meta$fullID <- NULL

#-------------------------------------------------------------------------------
# Create otu_table
species_reads_long <- master_long_df_filtered %>%
  select(asv, reads, fullID) %>%
  mutate(across(c(asv, fullID), as.factor)) %>%
  distinct()

species_reads_wide <- species_reads_long %>% 
  pivot_wider(names_from = fullID, 
              values_from = reads, 
              values_fill = list(reads = 0)) %>% 
  as.data.frame()

#order names alphabetically
species_reads_wide<- species_reads_wide[order(species_reads_wide$asv),] #order rows alphabetically
species_reads_wide <- species_reads_wide[,order(colnames(species_reads_wide))] #order cols alphabetically

#update row names to ASVs
rownames(species_reads_wide) <- NULL
species_reads_wide <- as.data.frame(species_reads_wide)
rownames(species_reads_wide) <- species_reads_wide$asv
species_reads_wide$asv <- NULL

print(dim(species_reads_wide))

#-------------------------------------------------------------------------------
# create otu_table (taxa level) 
species_reads_long_taxa <- master_long_df_filtered %>%
  select(asv, taxa_name_final, reads, fullID) %>%
  mutate(across(c(asv, taxa_name_final, fullID), as.factor)) %>%
  distinct()

species_reads_long_taxa_grouped <- species_reads_long_taxa %>% 
  group_by(taxa_name_final, fullID) %>% 
  summarize(reads = sum(reads, na.rm = TRUE), .groups = "drop")

species_reads_wide_taxa <- species_reads_long_taxa_grouped %>% 
  pivot_wider(names_from = fullID, 
              values_from = reads, 
              values_fill = list(reads = 0)) %>% 
  as.data.frame()

#order names alphabetically
species_reads_wide_taxa<- species_reads_wide_taxa[order(species_reads_wide_taxa$taxa_name_final),] #order rows alphabetically
species_reads_wide_taxa <- species_reads_wide_taxa[,order(colnames(species_reads_wide_taxa))] #order cols alphabetically

#update row names to ASVs
rownames(species_reads_wide_taxa) <- NULL
species_reads_wide_taxa <- as.data.frame(species_reads_wide_taxa)
rownames(species_reads_wide_taxa) <- species_reads_wide_taxa$taxa_name_final
species_reads_wide_taxa$taxa_name_final <- NULL
dim(species_reads_wide_taxa)

#-------------------------------------------------------------------------------
# build phyloseq
phylo_asv <- otu_table(species_reads_wide, taxa_are_rows=TRUE)
phylo_asv_taxa <- otu_table(species_reads_wide_taxa, taxa_are_rows=TRUE)
phylo_tax <- tax_table(taxa)
phylo_samples <- sample_data(meta)

# Check if colnames in asv object match meta rownames
list1<- colnames(phylo_asv_taxa)
list2 <- rownames(phylo_samples)
which(list1 != list2) #should be NULL

# Check if species names match in asv object and taxa object
list1<- rownames(phylo_asv_taxa)
list2 <- rownames(phylo_tax)
which(list1 != list2) #should be NULL

# final object
phylo_eDNA <- phyloseq(phylo_asv_taxa, phylo_tax, phylo_samples)

# Remove taxa wihout occurances from the data
phylo_eDNA <- subset_taxa(phylo_eDNA, taxa_sums(phylo_eDNA) > 0)

# save
saveRDS(phylo_eDNA, "Data/Processed/phylo_eDNA.RDS")