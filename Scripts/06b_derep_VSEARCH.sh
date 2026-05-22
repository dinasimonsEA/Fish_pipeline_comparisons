#!/bin/bash
set -euo pipefail

# merging, dereplicating, and chimera removal using VSEARCH de novo
# different order to dada2 - 1. merge → 2. filter → 3. dereplicate → 4.denoise → 5.chimera removal
# need to debug

# ==============================
# PATHS
# ==============================

base="/Workspace/Users/dina.simons@environment-agency.gov.uk/Fish_pipeline_comparisons/Data"

input_filtered="${base}/Temp/04_filterAndTrim"
manifest="${base}/Temp/04_filterAndTrim/04_manifest_filtered.csv"

work="${base}/Temp/06b_vsearch_denovo_outputs"
results="${base}/Results"

mkdir -p "$work" "$results"

# ==============================
# 1. MERGING
# ==============================

merge_dir="${work}/01_merged"
mkdir -p "$merge_dir"

echo "=== MERGING ==="

: > "$merge_dir/merged-reads.txt"
: > "$merge_dir/vsearch.merge.out"

tail -n +2 "$manifest" | while IFS=',' read -r idx sample fwd rev
do
    sample=$(echo "$sample" | sed 's/"//g')
    fwd=$(echo "$fwd" | sed 's/"//g')
    rev=$(echo "$rev" | sed 's/"//g')

    merged="${sample}_merged.fastq"

    vsearch \
        --fastq_mergepairs "${fwd}" \
        --reverse "${rev}" \
        --fastqout "${merge_dir}/${merged}" \
        --fastq_allowmergestagger \
        >> "${merge_dir}/vsearch.merge.out" 2>&1

    echo "$merged" >> "$merge_dir/merged-reads.txt"
    echo "$sample DONE"
done

echo "Merging complete."
echo "--------------------------------------"

# ==============================
# 2. FILTER
# ==============================

filter_dir="${work}/02_filtered"
mkdir -p "$filter_dir"

echo "=== FILTERING ==="

: > "$filter_dir/filtered-reads.txt"
: > "$filter_dir/filter.out"

while read merged
do
    base=$(basename "$merged" .fastq)
    filtered="${base}_filtered.fasta"

    vsearch \
        --fastq_filter "${merge_dir}/${merged}" \
        --fastq_maxee 2.0 \
        --relabel Filt \
        --fastaout "${filter_dir}/${filtered}" \
        >> "${filter_dir}/filter.out" 2>&1

    echo "$filtered" >> "$filter_dir/filtered-reads.txt"
    echo "$filtered DONE"
done < "${merge_dir}/merged-reads.txt"

# ==============================
# 3. DEREPLICATION
# ==============================

derep_dir="${work}/03_derep"
mkdir -p "$derep_dir"

echo "=== DEREPLICATION ==="

: > "$derep_dir/uniques-reads.txt"
: > "$derep_dir/uniques.out"

while read filtered
do
    base=$(basename "$filtered" .fasta)
    unique="${base}_uniques.fasta"

    vsearch \
        --derep_fulllength "${filter_dir}/${filtered}" \
        --sizeout \
        --relabel Uniq \
        --output "${derep_dir}/${unique}" \
        >> "${derep_dir}/uniques.out" 2>&1

    echo "$unique" >> "$derep_dir/uniques-reads.txt"
    echo "$unique DONE"
done < "${filter_dir}/filtered-reads.txt"

# ==============================
# 4. DENOISING (UNOISE)
# ==============================

denoise_dir="${work}/04_denoised"
mkdir -p "$denoise_dir"

echo "=== DENOISING ==="

: > "$denoise_dir/denoised-reads.txt"
: > "$denoise_dir/denoised.out"

while read unique
do
    base=$(basename "$unique" .fasta)
    denoised="${base}_denoised.fasta"

    vsearch \
        --cluster_unoise "${derep_dir}/${unique}" \
        --sizein \
        --sizeout \
        --centroids "${denoise_dir}/${denoised}" \
        >> "${denoise_dir}/denoised.out" 2>&1

    echo "$denoised" >> "$denoise_dir/denoised-reads.txt"
    echo "$denoised DONE"
done < "${derep_dir}/uniques-reads.txt"

# ==============================
# 5. CHIMERA REMOVAL
# ==============================

nochim_dir="${work}/05_nochim"
mkdir -p "$nochim_dir"

echo "=== CHIMERA REMOVAL ==="

: > "$nochim_dir/nochim-reads.txt"
: > "$nochim_dir/nochim.out"

while read denoised
do
    base=$(basename "$denoised" .fasta)
    nochim="${base}_nochim.fasta"

    vsearch \
        --uchime3_denovo "${denoise_dir}/${denoised}" \
        --sizein \
        --nonchimeras "${nochim_dir}/${nochim}" \
        --uchimeout "${nochim_dir}/${base}_uchime.txt" \
        >> "${nochim_dir}/nochim.out" 2>&1

    echo "$nochim" >> "$nochim_dir/nochim-reads.txt"
    echo "$nochim DONE"
done < "${denoise_dir}/denoised-reads.txt"

echo ="PIPELINE DONE"
