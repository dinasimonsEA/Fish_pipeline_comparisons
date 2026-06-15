# script to take tidied long data and turn into phyloseq object

# ----------------------------------------------------------------
# create tax_table()

# =========================
# STEP 1: MERGE TAXONOMY
# =========================

taxa_clean <- master_long_worms_df %>%
  mutate(
    kingdom = if_else(!is.na(AphiaID_worms), kingdom_worms, kingdom),
    phylum  = if_else(!is.na(AphiaID_worms), phylum_worms, phylum),
    class   = if_else(!is.na(AphiaID_worms), class_worms, class),
    order   = if_else(!is.na(AphiaID_worms), order_worms, order),
    family  = if_else(!is.na(AphiaID_worms), family_worms, family),
    genus   = if_else(!is.na(AphiaID_worms), genus_worms, genus),
    species = if_else(!is.na(AphiaID_worms), species_worms, species),
    aphiaid = AphiaID_worms
  )

# =========================
# STEP 2: STRIP ALL Rle STRUCTURES (CRITICAL)
# =========================

taxa <- taxa_clean %>%
  select(
    aphiaid,
    taxa_name_clean,
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

# convert back to tibble
taxa <- as_tibble(taxa)

# =========================
# STEP 3: PRIORITY LOGIC
# =========================

taxa <- taxa %>%
  mutate(
    priority = case_when(
      !is.na(aphiaid) & !is.na(kingdom) ~ 3,
      !is.na(kingdom) ~ 2,
      TRUE ~ 1
    )
  )

# =========================
# STEP 4: REMOVE DUPLICATES (SAFE)
# =========================

taxa <- taxa %>%
  group_by(taxa_name_clean) %>%
  arrange(desc(priority), .by_group = TRUE) %>%
  slice_head(n = 1) %>%   # safer than slice(1)
  ungroup() %>%
  select(-priority)

# =========================
# STEP 5: FINAL CLEANUP
# =========================

taxa <- taxa %>%
  distinct() %>%
  as_tibble()

rownames(taxa) <- NULL

# =========================
# STEP 6: SORT FOR PHYLOSEQ
# =========================

taxa$taxa_name_clean <- as.factor(taxa$taxa_name_clean)
taxa <- taxa[order(taxa$taxa_name_clean), ]

# =========================
# STEP 7: TAX MATRIX
# =========================

tax_matrix <- taxa %>%
  column_to_rownames("taxa_name_clean") %>%
  as.matrix()

# =========================
# CHECKS
# =========================

cat("Duplicate taxa:", any(duplicated(rownames(tax_matrix))), "\n")
cat("Missing kingdom:", sum(is.na(taxa$kingdom)), "\n")

# turn into phyloseq format
taxa$taxa_name_clean <- as.factor(taxa$taxa_name_clean)
taxa <- taxa[order(taxa$taxa_name_clean),] #order alphabetically to match other ps elements
rownames(taxa) <- NULL # reset row names
rownames(taxa) <-taxa$taxa_name_clean
taxa <- as.matrix(taxa)

#-------------------------------------------------------------------------------
# create sample_data()
meta <- unique(subset(master_long_df, select = c("sampleid",
                                             "denoise_method",
                                             "taxonomy_method",
                                             "dataset",
                                             "fullID"
)))

meta <- meta[order(meta$fullID),] #order rows alphabetically
rownames(meta) <-meta$fullID

#-------------------------------------------------------------------------------
# Create otu_table
species_reads_long <- master_long_df %>%
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
species_reads_wide <- species_reads_wide %>% column_to_rownames(var="asv")

print(dim(species_reads_wide))

#-------------------------------------------------------------------------------
# create otu_table (taxa level) 
species_reads_long_taxa <- master_long_df %>%
  select(asv, taxa_name_clean, reads, fullID) %>%
  mutate(across(c(asv, taxa_name_clean, fullID), as.factor)) %>%
  distinct()

species_reads_long_taxa_grouped <- species_reads_long_taxa %>% 
  group_by(taxa_name_clean, fullID) %>% 
  summarize(reads = sum(reads, na.rm = TRUE), .groups = "drop")

species_reads_wide_taxa <- species_reads_long_taxa_grouped %>% 
  pivot_wider(names_from = fullID, 
              values_from = reads, 
              values_fill = list(reads = 0)) %>% 
  as.data.frame()

#order names alphabetically
species_reads_wide_taxa<- species_reads_wide_taxa[order(species_reads_wide_taxa$taxa_name_clean),] #order rows alphabetically
species_reads_wide_taxa <- species_reads_wide_taxa[,order(colnames(species_reads_wide_taxa))] #order cols alphabetically

#update row names to ASVs
rownames(species_reads_wide_taxa) <- NULL
species_reads_wide_taxa <- species_reads_wide_taxa %>% column_to_rownames(var="taxa_name_clean")

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