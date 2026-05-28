# DOES NOT RUN: script for condensing taxonomy for BLAST outputs using LCA
# part of old code only

#-------------------------------
#previous script
#Need to format the accessions database first though.

%sh
#Filter large accession db
# directly downloaded larger accession dataset and uploaded to lab zone
#filtering to reduce size. created in sql in next script.

test_data_name="Marchamley"

cat Data/Processed/${test_data_name}/08b_ASVs_*_blast_Metafishlib.txt | cut -f2 | sort | uniq > Data/Databases/${test_data_name}_accessions.txt

zgrep -Ff Data/Databases/${test_data_name}_accessions.txt /Volumes/prd_dash_lab/ea_csg_research_edna_restricted/shared_external_volume/Fish_pipeline_comparison/Databases/nucl_gb.accession2taxid.gz > Data/Databases/${test_data_name}_subset.txt

# Create sql file: Takes a text file of accession → taxid mappings and converts it into a fast, searchable SQLite database

# convert subset to SQLite database
read.accession2taxid(
  taxaFiles = file.path("Data", "Databases", paste0(test_data_name, "_subset.txt")),
  sqlFile   = file.path("Data", "Databases", paste0(test_data_name, "_small_accessionTaxa.sql")),
  overwrite = TRUE
)

# store path
sql <- file.path("Data", "Databases", paste0(test_data_name, "_small_accessionTaxa.sql"))

## more standard options but they don't work in DASH - other code in main notebook
#options(timeout = 6000)  # increase to ~100 minutes
#taxonomizr::prepareDatabase('accessionTaxa.sql')

#direct download
#wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz

#prepareDatabase(
  #sqlFile = "/Volumes/prd_dash_lab/ea_csg_research_edna_restricted/shared_external_volume/Fish_pipeline_comparison/Databases/accessionTaxa.sql",
  #accessionTaxaFiles = "/Volumes/prd_dash_lab/ea_csg_research_edna_restricted/shared_external_volume/Fish_pipeline_comparison/Databases/nucl_gb.accession2taxid.gz")

# Getting info ------------------------------------------
dada2_accessions <- BLAST_DADA2_output$sseqid 
dada2_taxaId<-accessionToTaxa(dada2_accessions, sql)
print("Preview of taxa IDs for DADA2:")
head(dada2_taxaId)

VSEARCH_accessions <- BLAST_VSEARCH_output$sseqid 
VSEARCH_taxaId<-accessionToTaxa(VSEARCH_accessions, sql)
print("Preview of taxa IDs for VSEARCH:")
head(VSEARCH_taxaId)

## get taxonomy for names provided from BLAST
data2_taxonomy <- getTaxonomy(taxaId,'accessionTaxa.sql')

#taxaId_NCBI<-accessionToTaxa(accessions_NCBI,"accessionTaxa.sql")
#print("Preview of taxa IDs for NCBI:")
#print(head(taxaId_NCBI))

### get the taxonomy for those IDs
taxa_DADA2 <- getTaxonomy(taxaId,'accessionTaxa.sql')
print("Preview of taxa for DADA2:")
print(head(taxa_DADA2))

## find common names 
common_DADA2 <- getCommon(taxaId,'accessionTaxa.sql')
print("Preview of common taxa names for DADA2:")
print(head(common_DADA2))

## Condense multiple taxonomic assignments to their most recent common branch
taxa_LCA_DADA2 <- condenseTaxa(taxa_DADA2)

#taxa_LCA_NCBI <- condenseTaxa(common_NCBI)

# join all taxonomy data
#-----------------------------------------------------------------------

# -----------------------------
# USER SETTINGS
# -----------------------------
processed_dir <- file.path(path, "Data", "Processed", test_data_name)

blast_file_DADA2 <- file.path(processed_dir, "08b_ASVs_DADA2_blast_Metafishlib.txt")
blast_file_vsearch <- file.path(processed_dir, "08b_ASVs_VSEARCH_blast_Metafishlib.txt")

accession_sql <- file.path("Data","Databases", paste0(test_data_name, "_small_accessionTaxa.sql"))

nodes_file <- "nodes.dmp"
names_file <- "names.dmp"   # optional

# filtering thresholds
top_fraction <- 0.95   # keep hits within 95% of best bitscore
min_pident   <- 90     # optional, set NULL to skip

# -----------------------------
# LOAD BLAST DATA
# -----------------------------
blast <- read.delim(blast_file, stringsAsFactors = FALSE)

# EXPECTED columns: qseqid, sseqid, bitscore, pident
# adjust names if needed
colnames(blast)[1:4] <- c("qseqid","sseqid","pident","bitscore")

# -----------------------------
# FILTER BLAST HITS
# -----------------------------
blast <- blast %>%
  group_by(qseqid) %>%
  filter(bitscore >= top_fraction * max(bitscore)) %>%
  ungroup()

if (!is.null(min_pident)) {
  blast <- blast %>% filter(pident >= min_pident)
}

# -----------------------------
# ACCESSION → TAXID
# -----------------------------
blast$taxid <- accessionToTaxa(blast$sseqid, accession_sql)

blast <- blast %>% filter(!is.na(taxid))

# -----------------------------
# LOAD TAXONOMY (nodes)
# -----------------------------
nodes <- read.delim(nodes_file, sep="|", header=FALSE, stringsAsFactors=FALSE)
nodes <- nodes[, c(1,2)]
colnames(nodes) <- c("taxid","parent")

nodes$taxid  <- as.integer(trimws(nodes$taxid))
nodes$parent <- as.integer(trimws(nodes$parent))

# Speed lookup
parent_lookup <- setNames(nodes$parent, nodes$taxid)

# -----------------------------
# FUNCTIONS
# -----------------------------
get_lineage <- function(taxid, parent_lookup) {
  lineage <- c(taxid)
  
  while (!is.na(parent_lookup[as.character(taxid)]) &&
         parent_lookup[as.character(taxid)] != taxid) {
    
    taxid <- parent_lookup[as.character(taxid)]
    lineage <- c(lineage, taxid)
  }
  
  return(lineage)
}

get_lca <- function(taxids, parent_lookup) {
  lineages <- lapply(taxids, get_lineage, parent_lookup = parent_lookup)
  
  common <- Reduce(intersect, lineages)
  
  if (length(common) == 0) return(NA)
  
  # return lowest (closest to leaves)
  for (t in lineages[[1]]) {
    if (t %in% common) return(t)
  }
  
  return(NA)
}

# -----------------------------
# CALCULATE LCA PER ASV
# -----------------------------
lca_results <- blast %>%
  group_by(qseqid) %>%
  summarise(
    lca_taxid = get_lca(unique(taxid), parent_lookup)
  ) %>%
  ungroup()

# -----------------------------
# OPTIONAL: ADD TAXON NAMES
# -----------------------------
if (file.exists(names_file)) {
  
  names <- read.delim(names_file, sep="|", header=FALSE, stringsAsFactors=FALSE)
  names <- names[, c(1,2,4)]
  colnames(names) <- c("taxid","name","type")
  
  names <- names[names$type == "scientific name", ]
  
  names$taxid <- as.integer(trimws(names$taxid))
  names$name  <- trimws(names$name)
  
  lca_results <- merge(
    lca_results,
    names[, c("taxid","name")],
    by.x = "lca_taxid",
    by.y = "taxid",
    all.x = TRUE
  )
}

# -----------------------------
# SAVE OUTPUT
# -----------------------------
output_file <- file.path("Data","Processed", test_data_name, "LCA_results.txt")

write.table(
  lca_results,
  file = output_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("✅ LCA calculation complete\n")