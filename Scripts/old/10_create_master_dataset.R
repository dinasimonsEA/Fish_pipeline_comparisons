# script to create master datasets from the many files

# Import and merge ASV data ----------------------------------------

cat("\nImporting and merging ASVs...\n")

master_list <- list()

for (method in denoise_method) {
  for (ds in datasets) {
    
    cat("\nProcessing:", ds, "|", method, "\n")
    
    ds_name <- tolower(gsub("_.*", "", ds))
    method_name <- tolower(method)
    
    count_file <- file.path(base_path, ds, paste0("06_ASV_counts_", method, ".tsv"))
    fasta_file <- file.path(base_path, ds, paste0("06_ASV_seqs_", method, ".fasta"))
    
    # CHECK FILES EXIST
    if (!file.exists(count_file) | !file.exists(fasta_file)) {
      warning(paste("Missing files for", ds, method))
      next
    }
    
    # READ COUNT TABLE
    asv_count <- read.table(count_file, header = TRUE, row.names = 1, check.names = FALSE)
    
    ## Extract sequences from rownames
    asv_count$sequence <- toupper(trimws(rownames(asv_count)))
    
    # READ FASTA
    asv_fasta <- read.fasta(fasta_file)
    
    asv_fasta_df <- data.frame(
      ASV = names(asv_fasta),
      sequence = toupper(trimws(unlist(getSequence(asv_fasta, as.string = TRUE))))
    )
    
    # RESHAPE (WIDE → LONG)
    count_long <- stack(asv_count[, -ncol(asv_count)])  # exclude sequence column
    
    count_long$sequence <- rep(
      asv_count$sequence,
      times = ncol(asv_count) - 1
    )
    
    colnames(count_long) <- c("reads", "sampleID", "sequence")
    
    count_long <- count_long[, c("sequence", "sampleID", "reads")]
    
    ## Remove zero reads
    count_long <- count_long[count_long$reads > 0, ]
    
    # DEBUG MATCH CHECK
    match_count <- sum(count_long$sequence %in% asv_fasta_df$sequence)
    cat("Matching sequences:", match_count, "\n")
    
    if (match_count == 0) {
      warning(paste("No sequence matches for", ds, method))
      next
    }
    
    # MERGE (sequence → ASV)
    merged <- merge(count_long, asv_fasta_df, by = "sequence")
    
    if (nrow(merged) == 0) {
      warning(paste("Merge failed for", ds, method))
      next
    }
    
    # ADD METADATA
    merged$denoise_method <- method_name
    merged$dataset <- ds_name
    
    # FINAL FORMAT
    merged <- merged[, c("ASV", "sampleID", "reads", "sequence", "denoise_method", "dataset")]
    
    ## Store
    master_list[[length(master_list) + 1]] <- merged
  }
}

# COMBINE ALL DATASETS
if (length(master_list) > 0) {
  ASV_merged_dataset <- do.call(rbind, master_list)
  ASV_merged_dataset[c("ASV","sampleID","denoise_method","dataset")] <-
  lapply(ASV_merged_dataset[c("ASV","sampleID","denoise_method","dataset")], as.factor)
} else {
  stop("No data was successfully processed.")
}

# get taxonomy master ---------------------------------------

taxonomy_master_list <- list()

# -------------------------
# Helper: taxonomy string
# -------------------------
make_taxonomy_string <- function(df) {
  apply(df[, c("kingdom","phylum","class","order","family","genus","species")],
        1,
        function(x) paste(na.omit(x), collapse = ";"))
}

# -------------------------
# Flexible parser
# -------------------------
standardise_taxonomy <- function(df) {
  
  cn <- tolower(colnames(df))
  
  find_col <- function(pattern) {
    idx <- grep(pattern, cn)
    if (length(idx) > 0) idx[1] else NA
  }
  
  get_col <- function(idx) {
    if (!is.na(idx)) df[[idx]] else NA
  }
  
  # find ASV column robustly
  asv_idx <- grep("^asv$|asv_id", cn)
  if (length(asv_idx) == 0) asv_idx <- 1
  
  data.frame(
    ASV = df[[asv_idx[1]]],
    kingdom = get_col(find_col("kingdom")),
    phylum  = get_col(find_col("phylum")),
    class   = get_col(find_col("class")),
    order   = get_col(find_col("order")),
    family  = get_col(find_col("family")),
    genus   = get_col(find_col("genus")),
    species = get_col(find_col("species"))
  )
}

# -------------------------
# LCA parser
# -------------------------

standardise_lca <- function(df) {
  
  # initialise empty columns
  out <- data.frame(
    ASV = df$ASV_ID,
    kingdom = NA,
    phylum  = df$Phylum,
    class   = df$Class,
    order   = df$Order,
    family  = df$Family,
    genus   = df$Genus,
    species = df$Species,
    lca_rank = df$lca_rank,
    lca_name = df$lca_name
  )
  
  # fill missing taxonomy using LCA result
  for (i in seq_len(nrow(out))) {
    
    rank <- tolower(out$lca_rank[i])
    name <- out$lca_name[i]
    
    if (!is.na(rank) && !is.na(name)) {
      
      if (rank == "species") out$species[i] <- name
      if (rank == "genus")   out$genus[i]   <- name
      if (rank == "family")  out$family[i]  <- name
      if (rank == "order")   out$order[i]   <- name
      if (rank == "class")   out$class[i]   <- name
      if (rank == "phylum")  out$phylum[i]  <- name
    }
  }
  
  return(out[, c("ASV","kingdom","phylum","class","order","family","genus","species")])
}

# -------------------------
# MAPSEQ column names
# -------------------------
mapseq_cols <- c(
"ASV","hit_id","alignment_len","identity","aln_len2","mismatch","gapopen",
"qstart","qend","sstart","send","strand","kingdom","k_score","k_conf",
"phylum","p_score","p_conf","class","c_score","c_conf",
"order","o_score","o_conf","family","f_score","f_conf",
"genus","g_score","g_conf","species","s_score","s_conf"
)

# -------------------------
# MAIN LOOP
# -------------------------
for (ds in datasets) {
  for (method in denoise_method) {
    
    dataset_name <- tolower(ds)
    method_name  <- tolower(method)
    folder_path  <- file.path(base_path, ds)
    
    cat("\n===============================")
    cat("\nProcessing:", dataset_name, "|", method_name, "\n")
    
    # =========================
    # MAPSEQ
    # =========================
    mapseq_all <- list.files(folder_path,
      pattern = "mapseq.*riaz.*\\.mseq$",
      ignore.case = TRUE,
      full.names = TRUE
    )
    
    mapseq_file <- mapseq_all[
      grepl(paste0("_", method_name, "_"), basename(mapseq_all), ignore.case = TRUE)
    ]
    
    cat("\nMAPSEQ matches:", length(mapseq_file), "\n")
    print(mapseq_file)
    
    if (length(mapseq_file) == 1) {
      
      mapseq_data <- read.table(mapseq_file, header = FALSE)
      
      if (ncol(mapseq_data) == length(mapseq_cols)) {
        
        colnames(mapseq_data) <- mapseq_cols
        
        mapseq_tax <- mapseq_data[, c(
          "ASV","kingdom","phylum","class","order","family","genus","species"
        )]
        
        mapseq_tax$taxonomy <- make_taxonomy_string(mapseq_tax)
        mapseq_tax$taxonomy_method <- "mapseq"
        mapseq_tax$dataset <- dataset_name
        mapseq_tax$denoise_method <- method_name
        
        taxonomy_master_list[[length(taxonomy_master_list)+1]] <- mapseq_tax
      }
      
    } else warning(paste("MAPSEQ issue:", ds, method_name))

# =========================
# RDP
# =========================

standardise_rdp <- function(df) {
  
  # ✅ FIX: extract ASV from rownames
  ASV_ids <- rownames(df)
  
  cn <- tolower(colnames(df))
  
  # Detect taxonomy column
  tax_idx <- grep("tax", cn)
  
  if (length(tax_idx) > 0) {
    
    taxonomy_col <- df[[tax_idx[1]]]
    
    # Split taxonomy string
    split_tax <- strsplit(as.character(taxonomy_col), ";")
    
    extract_rank <- function(x, prefix) {
      val <- x[grep(prefix, x)]
      if (length(val) > 0) {
        return(sub(prefix, "", val[1]))
      } else {
        return(NA)
      }
    }
    
    out <- data.frame(
      ASV = ASV_ids,   # ✅ NOW CORRECT
      kingdom = sapply(split_tax, extract_rank, "k__"),
      phylum  = sapply(split_tax, extract_rank, "p__"),
      class   = sapply(split_tax, extract_rank, "c__"),
      order   = sapply(split_tax, extract_rank, "o__"),
      family  = sapply(split_tax, extract_rank, "f__"),
      genus   = sapply(split_tax, extract_rank, "g__"),
      species = sapply(split_tax, extract_rank, "s__")
    )
    
    return(out)
    
  } else {
    
    # fallback: structured columns, still use rownames
    out <- standardise_taxonomy(df)
    out$ASV <- ASV_ids   # ✅ overwrite with correct IDs
    return(out)
  }
}
   
# RDP (STRICT MATCH)
rdp_all <- list.files(folder_path,
  pattern = "assigned_taxonomy.*riaz.*\\.csv$",
  ignore.case = TRUE,
  full.names = TRUE
)

rdp_file <- rdp_all[
  grepl(paste0("_", method_name, "_"), basename(rdp_all), ignore.case = TRUE)
]

cat("\nRDP matches:", length(rdp_file), "\n")
print(rdp_file)

if (length(rdp_file) == 1) {
  
  rdp_data <- read.csv(rdp_file, row.names=1, check.names=FALSE)
  
  # USE FIXED PARSER
  rdp_tax <- standardise_rdp(rdp_data)
  
  cat("RDP ASVs:", length(unique(rdp_tax$ASV)), "\n")
  
  rdp_tax$taxonomy <- make_taxonomy_string(rdp_tax)
  rdp_tax$taxonomy_method <- "RDP"
  rdp_tax$dataset <- dataset_name
  rdp_tax$denoise_method <- method_name
  
  taxonomy_master_list[[length(taxonomy_master_list)+1]] <- rdp_tax
  
} else warning(paste("RDP problem:", ds, method_name))  
    
    # =========================
    # LCA
    # =========================
    lca_all <- list.files(folder_path,
      pattern = "LCA.*riaz.*\\.tsv$",
      ignore.case = TRUE,
      full.names = TRUE
    )
    
    lca_file <- lca_all[
      grepl(toupper(method_name), basename(lca_all))
    ]
    
    cat("\nLCA matches:", length(lca_file), "\n")
    print(lca_file)
    
    if (length(lca_file) == 1) {
      
      lca_data <- read.table(lca_file, header = TRUE, sep = "\t")
      
      # use fixed parser
      lca_tax <- standardise_lca(lca_data)
    
      lca_tax$taxonomy <- make_taxonomy_string(lca_tax)
      lca_tax$taxonomy_method <- "BLAST"
      lca_tax$dataset <- dataset_name
      lca_tax$denoise_method <- method_name
      
      taxonomy_master_list[[length(taxonomy_master_list)+1]] <- lca_tax
      
    } else warning(paste("LCA issue:", ds, method_name))
  }
}

# -------------------------
# COMBINE
# -------------------------
cat("===============================")
cat("\nCombining taxonomy...\n")

taxonomy_master_dataset <- do.call(rbind, taxonomy_master_list)

cols <- c("ASV","kingdom","phylum","class","order","family","genus","species",
          "taxonomy_method","dataset","denoise_method")

taxonomy_master_dataset[cols] <- lapply(taxonomy_master_dataset[cols], as.factor)

cat("Total taxonomy rows:", nrow(taxonomy_master_dataset), "\n")


# -------------------------
# REMOVE DUPLICATES
# -------------------------
before <- nrow(taxonomy_master_dataset)

taxonomy_master_dataset <- taxonomy_master_dataset[
  !duplicated(taxonomy_master_dataset[, c("ASV","taxonomy_method","dataset","denoise_method")]),
]

cat("Duplicates removed:", before - nrow(taxonomy_master_dataset), "\n")

# -------------------------
# MERGE WITH ASV DATA
# -------------------------
cat("\nMerging with ASV dataset...\n")

final_master_dataset <- merge(
  ASV_merged_dataset,
  taxonomy_master_dataset,
  by = c("ASV", "dataset", "denoise_method"),
  all.x = TRUE
)

cat("Final rows:", nrow(final_master_dataset), "\n")
cat("Missing taxonomy:", sum(is.na(final_master_dataset$taxonomy)), "\n")

# save
write.csv(final_master_dataset, "Data/Processed/master_long_data.csv")

