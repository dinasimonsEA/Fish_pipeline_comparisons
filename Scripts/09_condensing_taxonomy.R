# script for condensing taxonomy for BLAST outputs using LCA
# uses taxonomiser

# accession database-----------------------------------------

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
VSEARCH_accessions <- BLAST_VSEARCH_output$sseqid 

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

