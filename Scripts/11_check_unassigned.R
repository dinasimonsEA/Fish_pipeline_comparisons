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
  input_data$unassigned <- grepl(
    paste(unassigned_phrases, collapse = "|"),
    input_data[[level]],
    ignore.case = TRUE
  ) | is.na(input_data[[level]])

  # -------------------------------
  # ASV normalised by reads
  # -------------------------------
  
read_summary <- aggregate(
    reads ~ denoise_method + taxonomy_method + dataset,
    data = input_data,
    FUN = sum
  )

  unassigned_reads <- aggregate(
    reads ~ denoise_method + taxonomy_method + dataset,
    data = subset(input_data, unassigned),
    FUN = sum
  )

  read_df <- merge(
    read_summary,
    unassigned_reads,
    by = c("denoise_method", "taxonomy_method", "dataset"),
    all.x = TRUE,
    suffixes = c("_total", "_unassigned")
  )

  #read_df$unassigned[is.na(read_df$unassigned)] <- 0

  read_df <- data.frame(
    level = level,
    metric = "Reads",
    denoise_method = read_df$denoise_method,
    taxonomy_method = read_df$taxonomy_method,
    dataset = read_df$dataset,
    total = read_df$reads_total,
    unassigned = read_df$reads_unassigned
  )

  read_df$percent_unassigned <- round(100 * read_df$unassigned / read_df$total, 2)

  # -------------------------------
  # ASV LEVEL
  # -------------------------------
  asv_unique <- unique(
    input_data[, c("asv","denoise_method","taxonomy_method","dataset","unassigned")]
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
  # TAXA LEVEL
  # -------------------------------
  taxa_unique <- unique(
    input_data[, c(level,"denoise_method","taxonomy_method","dataset","unassigned")]
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
  summary_list[[level]] <- rbind(read_df, asv_df, taxa_df)
}

# -------------------------------
# FINAL TABLE
# -------------------------------
unassigned_summary_table <- do.call(rbind, summary_list)