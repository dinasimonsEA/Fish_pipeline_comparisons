# script to BLAST sequences (alignment-based)

# define paths
db_src="/Workspace/Users/dina.simons@environment-agency.gov.uk/Fish_pipeline_comparisons/Data/Databases/12S_fish_db"
db_tmp="/tmp/blastdb"
db_name="12S_fish_db"
query="/Workspace/Users/dina.simons@environment-agency.gov.uk/Fish_pipeline_comparisons/Data/Processed/06_ASV_seqs_DADA2.fasta"
output_MFL="/Workspace/Users/dina.simons@environment-agency.gov.uk/Fish_pipeline_comparisons/Data/Processed/08b_ASVs_DADA2_blast_Metafishlib.txt"
output_NCBI="/Workspace/Users/dina.simons@environment-agency.gov.uk/Fish_pipeline_comparisons/Data/Processed/08b_ASVs_DADA2_blast_NCBI.txt"

# copy BLAST database to local disk (doesn't work on DASH without doing this)
mkdir -p $db_tmp
cp ${db_src}/${db_name}.* $db_tmp/

# blast against meta fish lib
blastn -query $query -task blastn -db ${db_tmp}/${db_name} -evalue 10 -outfmt "6 qseqid saccver stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore" -perc_identity 97 -out $output_MFL 

# blast against NCBI database - timeout issue when pulling directly
# need to upload a version of the NCBI database to lab storage and link

# using RNA ref seqs only here 

#split -l 200 $query query_chunk_

#for f in query_chunk_*; do
  #blastn -query $f \
    #-db refseq_rna \
    #-remote \
    #-task blastn \
    #-evalue 10 \
    #-outfmt 6 \
    #-perc_identity 97 \
    #-out ${f}.out
#done

#cat query_chunk_*.out > $output_NCBI

# filter for 97% only if not already set in blastn
#awk '($3>97)' Data/Processed/08b_ASVs_blast_Metafishlib.txt > Data/Processed/08b_ASVs_blast_Metafishlib_97.txt