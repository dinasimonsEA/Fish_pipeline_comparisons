#!/bin/bash
set -euo pipefail

# ==============================
# CONFIG
# ==============================

test_data_name="Marchamley"

base="/Workspace/Users/dina.simons@environment-agency.gov.uk/Fish_pipeline_comparisons/Data"

processed_dir="${base}/Processed/${test_data_name}"
db_src="${base}/Databases/12S_fish_db"

# BLAST temp directory (local disk)
db_tmp="/tmp/blastdb"
db_name="12S_fish_db"

# input/output files
query="${processed_dir}/06_ASV_seqs_DADA2.fasta"
output_MFL="${processed_dir}/08b_ASVs_DADA2_blast_Metafishlib.txt"
output_NCBI="${processed_dir}/08b_ASVs_DADA2_blast_NCBI.txt"

# ==============================
# CHECKS
# ==============================

if [[ ! -f "$query" ]]; then
    echo "ERROR: Query FASTA not found: $query"
    exit 1
fi

if ! command -v blastn &> /dev/null; then
    echo "ERROR: blastn is not available in PATH"
    exit 1
fi

# ==============================
# COPY DATABASE TO LOCAL DISK
# ==============================

echo "Copying BLAST database to local disk..."

mkdir -p "$db_tmp"
cp "${db_src}/${db_name}."* "$db_tmp/"

# check DB copied
if [[ ! -f "${db_tmp}/${db_name}.nhr" ]]; then
    echo "ERROR: BLAST database copy failed"
    exit 1
fi

# ==============================
# RUN BLAST (Meta Fish Lib)
# ==============================

echo "Running BLAST against Meta Fish Lib..."

blastn \
    -query "$query" \
    -task blastn \
    -db "${db_tmp}/${db_name}" \
    -evalue 10 \
    -outfmt "6 qseqid saccver stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
    -perc_identity 97 \
    -out "$output_MFL"

echo "BLAST (Meta Fish Lib) complete -> $output_MFL"

# ==============================
# OPTIONAL: NCBI REMOTE BLAST
# ==============================

# Uncomment if you re-enable this

# split query into chunks to avoid timeout
# split -l 200 "$query" query_chunk_

# for f in query_chunk_*; do
#     blastn \
#         -query "$f" \
#         -db refseq_rna \
#         -remote \
#         -task blastn \
#         -evalue 10 \
#         -outfmt 6 \
#         -perc_identity 97 \
#         -out "${f}.out"
# done

# cat query_chunk_*.out > "$output_NCBI"

# ==============================
# OPTIONAL FILTERING
# ==============================

# filter >97% identity (if needed)
# awk '($3 > 97)' "$output_MFL" > "${processed_dir}/08b_ASVs_DADA2_blast_Metafishlib_97.txt"

echo "DONE"