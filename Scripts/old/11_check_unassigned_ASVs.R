# script to summarise ASVs not assigned

# -------------------------------
# DEFINE UNASSIGNED TERMS
# -------------------------------
unassigned_phrases <- c(
  "uncultured","unclassified","ncbi","pseudo","not assigned",
  "sar","cellular","environmental","eukaryote","eukaryota","na"
)

# -------------------------------
# CREATE UNASSIGNED FLAG
# -------------------------------
master_long_df$unassigned <- grepl(
  paste(unassigned_phrases, collapse = "|"),
  master_long_df$species,
  ignore.case = TRUE
) | is.na(master_long_df$species)

cat("\nUnassigned flag created...\n")

# -------------------------------
# ---- ROW-LEVEL SUMMARY -------
# -------------------------------
cat("\n--- ROW LEVEL SUMMARY ---\n")

row_summary <- aggregate(
  unassigned ~ denoise_method + taxonomy_method + dataset,
  data = master_long_df,
  FUN = function(x) c(
    total = length(x),
    unassigned = sum(x),
    percent = round(100 * sum(x) / length(x), 2)
  )
)

# Clean output (important)
row_summary_clean <- data.frame(
  denoise_method = row_summary$denoise_method,
  taxonomy_method = row_summary$taxonomy_method,
  dataset = row_summary$dataset,
  total = row_summary$unassigned[, "total"],
  unassigned = row_summary$unassigned[, "unassigned"],
  percent_unassigned = row_summary$unassigned[, "percent"]
)

print(row_summary_clean)

# -------------------------------
# ---- ASV-LEVEL SUMMARY -------
# (unique ASVs per method)
# -------------------------------
cat("\n--- ASV LEVEL SUMMARY ---\n")

asv_unique <- unique(
  master_long_df[, c("asv","denoise_method","taxonomy_method","dataset", "unassigned")]
)

asv_summary <- aggregate(
  unassigned ~ denoise_method + taxonomy_method + dataset,
  data = asv_unique,
  FUN = function(x) c(
    total_ASVs = length(x),
    unassigned_ASVs = sum(x),
    percent = round(100 * sum(x) / length(x), 2)
  )
)

asv_summary_clean <- data.frame(
  denoise_method = asv_summary$denoise_method,
  taxonomy_method = asv_summary$taxonomy_method,
  dataset = asv_summary$dataset,
  total_ASVs = asv_summary$unassigned[, "total_ASVs"],
  unassigned_ASVs = asv_summary$unassigned[, "unassigned_ASVs"],
  percent_unassigned = asv_summary$unassigned[, "percent"]
)

print(asv_summary_clean)

# -------------------------------
# ---- TAXA (SPECIES) SUMMARY ---
# -------------------------------
cat("\n--- TAXA (SPECIES) LEVEL SUMMARY ---\n")

# get unique species per method
species_unique <- unique(
  master_long_df[, c("species","denoise_method","taxonomy_method","dataset","unassigned")]
)

# summarise
species_summary <- aggregate(
  unassigned ~ denoise_method + taxonomy_method + dataset,
  data = species_unique,
  FUN = function(x) c(
    total_species = length(x),
    unassigned_species = sum(x),
    percent = round(100 * sum(x) / length(x), 2)
  )
)

# clean output
species_summary_clean <- data.frame(
  denoise_method = species_summary$denoise_method,
  taxonomy_method = species_summary$taxonomy_method,
  dataset = species_summary$dataset,
  total_species = species_summary$unassigned[, "total_species"],
  unassigned_species = species_summary$unassigned[, "unassigned_species"],
  percent_unassigned = species_summary$unassigned[, "percent"]
)

print(species_summary_clean)

# -------------------------------
# ---- OPTIONAL QUICK TABLE -----
# -------------------------------
cat("\n--- QUICK CROSS TAB ---\n")

print(
  table(
    master_long_df$denoise_method,
    master_long_df$taxonomy_method,
    master_long_df$dataset,
    master_long_df$unassigned
  )
)

# -------------------------------
# ---- OPTIONAL CLEAN SUBSET ----
# -------------------------------
unassigned_data <- master_long_df[master_long_df$unassigned, ]

cat("\nTotal unassigned rows:", nrow(unassigned_data), "\n")