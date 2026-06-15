# -------------------------------
# SETUP
# -------------------------------
tax_levels <- c("species", "genus", "family")

unassigned_phrases <- c(
  "uncultured","unclassified","ncbi","pseudo","not assigned",
  "sar","cellular","environmental","eukaryote","eukaryota","na"
)

# storage
summary_list <- list()

# -------------------------------
# LOOP
# -------------------------------
for (level in tax_levels) {

  cat("Processing:", level, "\n")

  # -------------------------------
  # FLAG UNASSIGNED
  # -------------------------------
  master_long_df$unassigned <- grepl(
    paste(unassigned_phrases, collapse = "|"),
    master_long_df[[level]],
    ignore.case = TRUE
  ) | is.na(master_long_df[[level]])

  # -------------------------------
  # 1. ROW LEVEL
  # -------------------------------
  row_summary <- aggregate(
    unassigned ~ denoise_method + taxonomy_method + dataset,
    data = master_long_df,
    FUN = function(x) c(
      total = length(x),
      unassigned = sum(x)
    )
  )

  row_df <- data.frame(
    level = level,
    metric = "row",
    denoise_method = row_summary$denoise_method,
    taxonomy_method = row_summary$taxonomy_method,
    dataset = row_summary$dataset,
    total = row_summary$unassigned[, "total"],
    unassigned = row_summary$unassigned[, "unassigned"]
  )

  row_df$percent_unassigned <- round(100 * row_df$unassigned / row_df$total, 2)

  # -------------------------------
  # 2. ASV LEVEL
  # -------------------------------
  asv_unique <- unique(
    master_long_df[, c("asv","denoise_method","taxonomy_method","dataset","unassigned")]
  )

  asv_summary <- aggregate(
    unassigned ~ denoise_method + taxonomy_method + dataset,
    data = asv_unique,
    FUN = function(x) c(
      total = length(x),
      unassigned = sum(x)
    )
  )

  asv_df <- data.frame(
    level = level,
    metric = "ASV",
    denoise_method = asv_summary$denoise_method,
    taxonomy_method = asv_summary$taxonomy_method,
    dataset = asv_summary$dataset,
    total = asv_summary$unassigned[, "total"],
    unassigned = asv_summary$unassigned[, "unassigned"]
  )

  asv_df$percent_unassigned <- round(100 * asv_df$unassigned / asv_df$total, 2)

  # -------------------------------
  # 3. TAXA LEVEL
  # -------------------------------
  taxa_unique <- unique(
    master_long_df[, c(level,"denoise_method","taxonomy_method","dataset","unassigned")]
  )

  taxa_summary <- aggregate(
    unassigned ~ denoise_method + taxonomy_method + dataset,
    data = taxa_unique,
    FUN = function(x) c(
      total = length(x),
      unassigned = sum(x)
    )
  )

  taxa_df <- data.frame(
    level = level,
    metric = "taxa",
    denoise_method = taxa_summary$denoise_method,
    taxonomy_method = taxa_summary$taxonomy_method,
    dataset = taxa_summary$dataset,
    total = taxa_summary$unassigned[, "total"],
    unassigned = taxa_summary$unassigned[, "unassigned"]
  )

  taxa_df$percent_unassigned <- round(100 * taxa_df$unassigned / taxa_df$total, 2)

  # -------------------------------
  # COMBINE
  # -------------------------------
  summary_list[[level]] <- rbind(row_df, asv_df, taxa_df)
}

# -------------------------------
# FINAL TABLE
# -------------------------------
unassigned_summary_table <- do.call(rbind, summary_list)

# -------------------------------
# EXPORT
# -------------------------------
write.csv(unassigned_summary_table, "Results/unassigned_summary_table.csv", row.names = FALSE)

cat("\nSummary table created and saved as 'unassigned_summary_table.csv'\n")