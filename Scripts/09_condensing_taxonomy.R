# script for condensing taxonomy for BLAST outputs using LCA
# uses taxonomiser

# get Database
options(timeout = 6000)  # increase to ~100 minutes
prepareDatabase('accessionTaxa.sql')

# get BLAST results
blastResults_MFL <- read.table('Data/Processed/08b_ASVs_blast_Metafishlib.txt',header=FALSE,stringsAsFactors=FALSE)
blastResults_NCBI <- read.table('Data/Processed/08b_ASVs_blast_NCBI.txt',header=FALSE,stringsAsFactors=FALSE)

# Assigning taxonomy

## Producing accession numbers
accessions_MFL<-sapply(strsplit(blastResults_MFL[,2],'\\|'),'[',4) ## grab the 4th |-separated field from the reference name in the second column
accessions_NCBI<-sapply(strsplit(blastResults_NCBI[,2],'\\|'),'[',4)

## Finding taxonomy for NCBI accession numbers
### convert NCBI accession numbers to taxonomic IDs
taxaId_MFL<-accessionToTaxa(accessions_MFL,"accessionTaxa.sql")
print("Preview of taxa IDs for MFL:")
print(head(taxaId_MFL))

taxaId_NCBI<-accessionToTaxa(accessions_NCBI,"accessionTaxa.sql")
print("Preview of taxa IDs for NCBI:")
print(head(taxaId_NCBI))

### get the taxonomy for those IDs
taxa_MFL <- getTaxonomy(taxaId,'accessionTaxa.sql')
print("Preview of taxa for MFL:")
print(head(taxa_MFL))

taxa_NCBI <- getTaxonomy(taxaId,'accessionTaxa.sql')
print("Preview of taxa for NCBI:")
print(head(taxa_NCBI))

## find common names 
common_MFL <- getCommon(taxaId,'accessionTaxa.sql')
print("Preview of common taxa names for MFL:")
print(head(common_MFL))

common_NCBI <- getCommon(taxaId,'accessionTaxa.sql')
print("Preview of common taxa names for NCBI:")
print(head(common_NCBI))

## Condense multiple taxonomic assignments to their most recent common branch
taxa_LCA_MFL <- condenseTaxa(taxa_MFL)
taxa_LCA_NCBI <- condenseTaxa(common_NCBI)

# join all taxonomy data

