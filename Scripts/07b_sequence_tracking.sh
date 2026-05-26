#!/bin/bash
set -euo pipefail
shopt -s nullglob

# ==============================
# CONFIG
# ==============================

test_data_name="Marchamley"

base="/Workspace/Users/dina.simons@environment-agency.gov.uk/Fish_pipeline_comparisons/Data"
work="${base}/Temp/${test_data_name}/06b_vsearch_denovo_outputs"

mkdir -p "$work"

out="${work}/07b_vsearch_summary.tsv"

echo -e "sample\tmerged\tfiltered\tderep\tdenoised\tnochim" > "$out"

# ==============================
# FUNCTIONS
# ==============================

count_fastq() {
    [ -s "$1" ] && echo $(( $(wc -l < "$1") / 4 )) || echo 0
}

count_fasta() {
    [ -s "$1" ] && grep -c "^>" "$1" || echo 0
}

count_sizes() {
    if [ -s "$1" ]; then
        grep -o "size=[0-9]*" "$1" 2>/dev/null | cut -d= -f2 | awk '{s+=$1} END{print s+0}'
    else
        echo 0
    fi
}

echo "Building summary..."

# ==============================
# MAIN LOOP
# ==============================

for f in "${work}/01_merged/"*_merged.fastq
do
    # skip if no files found
    [ -e "$f" ] || continue

    sample=$(basename "$f" _merged.fastq)

    merged=$(count_fastq "$f")

    filt_file="${work}/02_filtered/${sample}_merged_filtered.fasta"
    derep_file="${work}/03_derep/${sample}_merged_filtered_uniques.fasta"
    denoise_file="${work}/04_denoised/${sample}_merged_filtered_uniques_denoised.fasta"
    nochim_file="${work}/05_nochim/${sample}_merged_filtered_uniques_denoised_nochim.fasta"

    filtered=0
    derep=0
    denoised=0
    nochim=0

    [ -f "$filt_file" ] && filtered=$(count_fasta "$filt_file")
    [ -f "$derep_file" ] && derep=$(count_sizes "$derep_file")
    [ -f "$denoise_file" ] && denoised=$(count_sizes "$denoise_file")
    [ -f "$nochim_file" ] && nochim=$(count_sizes "$nochim_file")

    echo -e "$sample\t$merged\t$filtered\t$derep\t$denoised\t$nochim" >> "$out"
done

echo "Done -> $out"

# ==============================
# PRINT SUMMARY
# ==============================

echo ""
echo "=== SUMMARY ==="
column -t "$out"