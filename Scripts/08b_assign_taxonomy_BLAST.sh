# script to BLAST sequences (alignment-based)

# to make BLAST formatted database from curated Meta Fish Lb
makeblastdb -in Data/Databases/Meta-fish-lib/references.12s.miya.qiime2.v268.fasta -dbtype nucl -out Data/Databases/12S_fish_db/12S_fish_db

# define paths
db_src="/Workspace/Users/dina.simons@environment-agency.gov.uk/Fish_pipeline_comparisons/Data/Databases/12S_fish_db"
db_tmp="/tmp/blastdb"
db_name="12S_fish_db"
query="/Workspace/Users/dina.simons@environment-agency.gov.uk/Fish_pipeline_comparisons/Data/Processed/06_ASV_seqs_DADA2.fasta"
output_MFL="/Workspace/Users/dina.simons@environment-agency.gov.uk/Fish_pipeline_comparisons/Data/Processed/08b_ASVs_blast_Metafishlib.txt"
output_NCBI="/Workspace/Users/dina.simons@environment-agency.gov.uk/Fish_pipeline_comparisons/Data/Processed/08b_ASVs_blast_NCBI.txt"

# copy BLAST database to local disk (doesn't work on DASH without doing this)
mkdir -p $db_tmp

cp ${db_src}/${db_name}.* $db_tmp/

# sanity check
ls -lh $db_tmp

# blast against meta fish lib
blastn -query $query -task blastn -db ${db_tmp}/${db_name} -evalue 10 -outfmt 6 -perc_identity 97 -out $output_MFL

# blast against NCBI database - usinf RNA ref seqs only here 
blastn -query $query -db refseq_rna -remote -task blastn -evalue 10 -outfmt 6 -perc_identity 97 -out $output_NCBI 

# filter for 97% only if not already set in blastn
#awk '($3>97)' Data/Processed/08b_ASVs_blast_Metafishlib.txt > Data/Processed/08b_ASVs_blast_Metafishlib_97.txt