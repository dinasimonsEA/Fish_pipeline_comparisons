# merging, dereplicating, and chimera removal using VSEARCH de novo
# different order to dada2 - 1. merge → 2. filter → 3. dereplicate → 4.denoise → 5.chimera removal

# 1. Merge pairs ---------------------

# Input directory
input_dir="/Workspace/Users/dina.simons@environment-agency.gov.uk/Fish_pipeline_comparisons/Data/Temp/filterAndTrim"

# Output directory
output_processed_dir="/Workspace/Users/dina.simons@environment-agency.gov.uk/Fish_pipeline_comparisons/Data/Processed"
output_temp_dir_merge="/Workspace/Users/dina.simons@environment-agency.gov.uk/Fish_pipeline_comparisons/Data/Temp/06b_vsearch_denovo_outputs/merged"

# make dirs
mkdir -p "$output_temp_dir_merge"
mkdir -p "$output_processed_dir"

# Identify manifest of filtered samples
manifest_filtered="/Workspace/Users/dina.simons@environment-agency.gov.uk/Fish_pipeline_comparisons/Data/Temp/manifest_filtered.csv"

# check files exist in input
ls "$input_dir"/*.fastq.gz || { echo "No FASTQ files found"; exit 1; }

# start
echo "Running vsearch merging"
echo "--------------------------------------"

# Clear previous outputs
rm -f "$output_temp_dir_merge"/vsearch.merge.out
rm -f "$output_temp_dir_merge"/merged-reads.txt

# loop
tail -n +2 "$manifest_filtered" | while IFS=',' read -r idx sample fwd rev
do
    # clean labels
    sample=$(echo "$sample" | sed 's/"//g')
    fwd=$(echo "$fwd" | sed 's/"//g')
    rev=$(echo "$rev" | sed 's/"//g')

    merged="${sample}_merged.fastq"

    echo "Merging: $sample"

    vsearch \
        --fastq_mergepairs "$fwd" \
        --reverse "$rev" \
        --fastqout "$output_temp_dir_merge/$merged" \
        --fastq_allowmergestagger \
        >> "$output_temp_dir_merge/vsearch.merge.out" 2>&1

    echo "${merged}" >> "$output_temp_dir_merge/merged-reads.txt"
    echo "${sample} DONE"
    echo "--------------------------------------"
done

# stats
grep -A 3 "Merging reads" "$output_temp_dir_merge/vsearch.merge.out" | grep "Pairs" > "$output_temp_dir_merge/pairs.column.txt"
grep -A 3 "Merging reads" "$output_temp_dir_merge/vsearch.merge.out" | grep "Merged" > "$output_temp_dir_merge/merged.column.txt"
grep -A 3 "Merging reads" "$output_temp_dir_merge/vsearch.merge.out" | grep "Not merged" > "$output_temp_dir_merge/not-merged.column.txt"

# summary table
paste "$output_temp_dir_merge/merged-reads.txt" \
      "$output_temp_dir_merge/pairs.column.txt" \
      "$output_temp_dir_merge/merged.column.txt" \
      "$output_temp_dir_merge/not-merged.column.txt" \
| awk '{print "merge", $0}' > "$output_temp_dir_merge/vsearch.merge.summary.txt"

# cleanup
rm "$output_temp_dir_merge"/*.column.txt

echo "Merging complete."
echo "--------------------------------------"

#2. Filter errors -------------------------------

# Input directory
input_dir="/Workspace/Users/dina.simons@environment-agency.gov.uk/Fish_pipeline_comparisons/Data/Temp/06b_vsearch_denovo_outputs/merged"

# Output directory
output_temp_dir_filtered="/Workspace/Users/dina.simons@environment-agency.gov.uk/Fish_pipeline_comparisons/Data/Temp/06b_vsearch_denovo_outputs/filtered"

mkdir -p "$output_temp_dir_filtered"

# start
echo "Running vsearch to filter errors"
echo "--------------------------------------"

for merged in `cat ${input_dir}/merged-reads.txt`
do
	filtered=${merged//merged.fastq/filtered.fasta}

	vsearch -fastq_filter ${input_dir}/${merged} -fastq_maxee 1.0 -relabel Filt -fastaout ${output_temp_dir_filtered}/${filtered} \
		>> ${output_temp_dir_filtered}/filter.out 2>&1

	fastq=$(echo ${filtered} | sed 's/.fasta.gz/.fasta/')
    
	mv ${output_temp_dir_filtered}/${filtered} ${output_temp_dir_filtered}/${fastq}
	echo ${fastq} >> ${output_temp_dir_filtered}/filtered-reads.txt
	echo ${fastq} "DONE"
    echo "--------------------------------------"
done

paste ${output_temp_dir_filtered}/filtered-reads.txt <(grep "sequences kept" ${output_temp_dir_filtered}/filter.out) \
	| awk '{print "filter", $0}' > ${output_temp_dir_filtered}/filter.summary.txt

# 3. Dereplicate -----------------------------------------------

# Input directory
input_dir="/Workspace/Users/dina.simons@environment-agency.gov.uk/Fish_pipeline_comparisons/Data/Temp/06b_vsearch_denovo_outputs/filtered"

# Output directory
output_temp_dir_derep="/Workspace/Users/dina.simons@environment-agency.gov.uk/Fish_pipeline_comparisons/Data/Temp/06b_vsearch_denovo_outputs/derep"

mkdir -p "$output_temp_dir_derep"

#start
echo "Running vsearch to dereplicate"
echo "--------------------------------------"

for filtered in `cat ${input_dir}/filtered-reads.txt`
do
	uniques=${filtered//filtered.fasta/uniques.fasta}
	vsearch -derep_fulllength ${input_dir}/${filtered} -sizeout -relabel Uniq \
		-output ${output_temp_dir_derep}/${uniques} >> ${output_temp_dir_derep}/uniques.out 2>&1
	fastq=$(echo ${uniques} | sed 's/.fasta.gz/.fasta/')

    # Rename .fasta.gz → .fasta
	mv ${output_temp_dir_derep}/${uniques} ${output_temp_dir_derep}/${fastq}

	echo ${fastq} >> ${output_temp_dir_derep}/uniques-reads.txt
	echo ${fastq} "DONE"
    echo "--------------------------------------"
done

paste ${output_temp_dir_derep}/uniques-reads.txt <(grep "unique sequences" ${output_temp_dir_derep}/uniques.out) \
	| awk '{print "derep", $0}' > ${output_temp_dir_derep}/uniques.summary.txt

# 4. Denoise -------------------------------------------------------

# Input directory
input_dir="/Workspace/Users/dina.simons@environment-agency.gov.uk/Fish_pipeline_comparisons/Data/Temp/06b_vsearch_denovo_outputs/derep"

# Output directory
output_temp_dir_denoised="/Workspace/Users/dina.simons@environment-agency.gov.uk/Fish_pipeline_comparisons/Data/Temp/06b_vsearch_denovo_outputs/denoised"

mkdir -p "$output_temp_dir_denoised"

#start
echo "Running vsearch to denoise"
echo "--------------------------------------"

for uniques in `cat ${input_dir}/uniques-reads.txt`
do
	# Define output filename
    denoised=${uniques//uniques.fasta/denoised.fasta}
	
    #Run denoise
    vsearch --cluster_unoise ${input_dir}/${uniques} --sizein --sizeout --centroids ${output_temp_dir_denoised}/${denoised}
            >> ${output_temp_dir_denoised}/denoised.out 2>&1
	
    # Rename .fasta.gz → .fasta
    fasta=$(echo ${denoised} | sed 's/.fasta.gz/.fasta/')
    mv ${output_temp_dir_denoised}/${denoised} ${output_temp_dir_denoised}/${fasta}

    # Track outputs
    echo ${fasta}.gz >> ${output_temp_dir_denoised}/denoised-reads.txt
    echo ${fasta}.gz "DONE"
    echo "--------------------------------------"

done

paste ${output_temp_dir_denoised}/denoised-reads.txt <(grep "unique sequences" ${output_temp_dir_denoised}/denoised.out) \
	| awk '{print "denoised", $0}' > ${output_temp_dir_denoised}/denoised.summary.txt


# 5. Chimera removal ------------------------------------------------

# Input directory
input_dir="/Workspace/Users/dina.simons@environment-agency.gov.uk/Fish_pipeline_comparisons/Data/Temp/06b_vsearch_denovo_outputs/denoised"

# Output directory
output_temp_dir_nochim="/Workspace/Users/dina.simons@environment-agency.gov.uk/Fish_pipeline_comparisons/Data/Temp/06b_vsearch_denovo_outputs/nochim"

mkdir -p "$output_temp_dir_nochim"

#start
echo "Running vsearch to remove chimeras"
echo "--------------------------------------"

for denoised in `cat ${input_dir}/denoised-reads.txt`
do
	# Define output filename
    nochim=${denoised//denoised.fasta/nochim.fasta}
    basename=$(basename ${denoised} .fasta.gz)
	
    #Run chimera removal
    vsearch \
        --uchime3_denovo ${input_dir}/${denoised} --sizein --nonchimeras ${output_temp_dir_nochim}/${nochim} --uchimeout ${output_temp_dir_nochim}/${basename}_uchime.txt \
        >> ${output_temp_dir_nochim}/nochim.out 2>&1
	
    # Rename .fasta.gz → .fasta
    fasta=$(echo ${nochim} | sed 's/.fasta.gz/.fasta/')
    mv ${output_temp_dir_nochim}/${nochim} ${output_temp_dir_nochim}/${fasta}

    # Track outputs
    echo ${fasta}.gz >> ${output_temp_dir_nochim}/nochim-reads.txt
    echo ${fasta}.gz "DONE"
    echo "--------------------------------------"

done

paste ${output_temp_dir_nochim}/nochim-reads.txt <(grep "unique sequences" ${output_temp_dir_nochim}/nochim.out) \
	| awk '{print "nochim", $0}' > ${output_temp_dir_nochim}/nochim.summary.txt


# join the summary outputs
cat /Workspace/Users/dina.simons@environment-agency.gov.uk/Fish_pipeline_comparisons/Data/Temp/06b_vsearch_denovo_outputs/*/*.summary.txt | sort -k1,1 > ${output_processed_dir}/06b_tracking_sequences_VSEARCH.txt
