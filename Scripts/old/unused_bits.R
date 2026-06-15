# DOES NOT RUN AS STONEALONE

#Import ASV data to individual objects -----------
for (method in denoise_method) {
  for (ds in datasets) {
    
    ds_name <- tolower(gsub("_.*", "", ds))
    method_name <- tolower(method)
    
    count_file <- file.path(base_path, ds, paste0("06_ASV_counts_", method, ".tsv"))
    fasta_file <- file.path(base_path, ds, paste0("06_ASV_seqs_", method, ".fasta"))
    
    asv_count <- read.table(count_file)
    asv_fasta <- read.fasta(fasta_file)
    
    asv_fasta_df <- data.frame(
      ASV = names(asv_fasta),
      Seqs = unlist(getSequence(asv_fasta, as.string = TRUE))
    )
    
    asv_fasta_df$Seqs <- toupper(asv_fasta_df$Seqs)
    
    assign(paste0(ds_name, "_asv_count_", method_name), asv_count)
    assign(paste0(ds_name, "_asv_fasta_", method_name), asv_fasta)
    assign(paste0(ds_name, "_asv_fasta_", method_name, "_df"), asv_fasta_df)
  }
}

# get taxonomy files individually ---------------

#blast_cols <- c(
  #"qseqid","saccver","stitle","pident","length",
  #"mismatch","gapopen","qstart","qend",
  #"sstart","send","evalue","bitscore"
#)

mapseq_cols <- c(
"ASV","hit_id","alignment_len","identity","aln_len2","mismatch","gapopen",
"qstart","qend","sstart","send","strand","kingdom","k_score","k_conf",
"phylum","p_score","p_conf","class","c_score","c_conf",
"order","o_score","o_conf","family","f_score","f_conf",
"genus","g_score","g_conf","species","s_score","s_conf"
)

for (ds in datasets) {
  
  ds_name <- tolower(gsub("_.*", "", ds))
  folder_path <- file.path(base_path, ds)
  
  # BLAST
  #blast_file <- list.files(
    #folder_path,
    #pattern = "blast.*EA_riaz.*\\.txt$",
    #full.names = TRUE
  #)
  
  #if (length(blast_file) > 0) {
    #blast_data <- read.table(blast_file[1], header = FALSE, sep = "\t")
    #colnames(blast_data) <- blast_cols
    #assign(paste0(ds_name, "_blast_ea_riaz"), blast_data)
  #}
  
  # ASSIGNED TAXONOMY
  assigned_file <- list.files(
    folder_path,
    pattern = "assigned_taxonomy.*EA.*riaz.*\\.csv$",
    full.names = TRUE
  )
  
  if (length(assigned_file) > 0) {
    assigned_data <- read.csv(assigned_file[1])
    assign(paste0(ds_name, "_RDP_ea_riaz"), assigned_data)
  }
  
  # MAPSEQ
  mapseq_file <- list.files(
    folder_path,
    pattern = "mapseq.*EA_riaz.*\\.mseq$",
    full.names = TRUE
  )
  
if (length(mapseq_file) > 0) {
    mapseq_data <- read.table(mapseq_file[1], header = FALSE)
    
    # Only assign names if column count matches
    if (ncol(mapseq_data) == length(mapseq_cols)) {
      colnames(mapseq_data) <- mapseq_cols
    } else {
      warning(paste("Column mismatch in MAPSEQ file:", ds))
    }
    
    assign(paste0(ds_name, "_mapseq_ea_riaz"), mapseq_data)
  }
  
  # BLAST LCA / CONDENSED TAXONOMY
  lca_file <- list.files(
    folder_path,
    pattern = "LCA_results.*\\.tsv$",
    full.names = TRUE
  )
  
  if (length(lca_file) > 0) {
    lca_data <- read.table(lca_file[1], header = TRUE, sep = "\t")
    assign(paste0(ds_name, "_lca_taxonomy"), lca_data)
  }
}

# tidying functions

# generic paser 
def standardise_taxonomy(df):
    cols = {c.lower(): c for c in df.columns}

    def get(name):
        for k,v in cols.items():
            if name in k:
                return df[v]
        return pd.Series([None]*len(df))

    asv_col = next((v for k,v in cols.items() if "asv" in k), df.columns[0])

    return pd.DataFrame({
        "ASV": df[asv_col],
        "kingdom": get("kingdom"),
        "phylum": get("phylum"),
        "class": get("class"),
        "order": get("order"),
        "family": get("family"),
        "genus": get("genus"),
        "species": get("species")
    })

# LCA parser 
def standardise_lca(df):
    out = pd.DataFrame({
        "ASV": df["ASV_ID"],
        "kingdom": None,
        "phylum": df["Phylum"],
        "class": df["Class"],
        "order": df["Order"],
        "family": df["Family"],
        "genus": df["Genus"],
        "species": df["Species"],
        "lca_rank": df["lca_rank"],
        "lca_name": df["lca_name"]
    })

    # vectorised fill
    mask = out["lca_name"].notna() & out["lca_rank"].notna()

    for rank in ["species","genus","family","order","class","phylum"]:
        idx = (out["lca_rank"].str.lower() == rank) & mask
        out.loc[idx, rank] = out.loc[idx, "lca_name"]

    return out[["ASV","kingdom","phylum","class","order","family","genus","species"]]

# RDP parser
def standardise_rdp(df):
    
    df = df.copy()
    df["ASV"] = df.index

    tax_col = next((c for c in df.columns if "tax" in c.lower()), None)

    if tax_col:
        split_tax = df[tax_col].astype(str).str.split(";")

        def extract(prefix):
            return split_tax.apply(
                lambda x: next((i.replace(prefix,"") for i in x if prefix in i), None)
            )

        return pd.DataFrame({
            "sequence": df["sequence"],
            "kingdom": extract("k__"),
            "phylum": extract("p__"),
            "class": extract("c__"),
            "order": extract("o__"),
            "family": extract("f__"),
            "genus": extract("g__"),
            "species": extract("s__")
        })

    else:
        return standardise_taxonomy(df)