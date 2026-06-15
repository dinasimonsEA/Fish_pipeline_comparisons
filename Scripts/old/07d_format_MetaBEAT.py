import pandas as pd
from collections import OrderedDict
import re
import os

# =========================================================
# DATASETS TO PROCESS
# =========================================================

datasets = [
    "Windermere_2017",
    "Marchamley",
    "RingTrial_Sean"
]

# =========================================================
# LOOP THROUGH DATASETS
# =========================================================

for ds in datasets:

    print("\n" + "="*60)
    print(f"Processing dataset: {ds}")
    print("="*60)

    # =====================================================
    # FILE PATHS
    # =====================================================

    base_path = f"/Workspace/Users/dina.simons@environment-agency.gov.uk/Fish_pipeline_comparisons/Data/Processed/{ds}"

    count_file = f"{base_path}/06_ASV_counts_MetaBEAT_untidied.tsv"
    fasta_file = f"{base_path}/06_ASV_seqs_MetaBEAT_untidied.fasta"

    output_file = f"{base_path}/06_ASV_counts_MetaBEAT.tsv"
    clean_fasta_output = f"{base_path}/06_ASV_seqs_MetaBEAT.fasta"

    # =====================================================
    # CHECK FILES EXIST
    # =====================================================

    missing_files = []

    if not os.path.exists(count_file):
        missing_files.append(count_file)

    if not os.path.exists(fasta_file):
        missing_files.append(fasta_file)

    if len(missing_files) > 0:

        print("\nERROR: Missing input file(s):")

        for mf in missing_files:
            print(f" - {mf}")

        print(f"Skipping dataset: {ds}")
        continue

    # =====================================================
    # READ FASTA FILE
    # Create mapping:
    # original ASV_ID -> sequence
    # =====================================================

    seq_dict = OrderedDict()

    try:

        with open(fasta_file, "r") as f:

            current_id = None
            current_seq = []

            for line in f:

                line = line.strip()

                if line.startswith(">"):

                    # Save previous sequence
                    if current_id is not None:
                        seq_dict[current_id] = "".join(current_seq)

                    # Start new sequence
                    current_id = line[1:]
                    current_seq = []

                else:
                    current_seq.append(line)

            # Save final sequence
            if current_id is not None:
                seq_dict[current_id] = "".join(current_seq)

        print(f"Loaded {len(seq_dict)} sequences from FASTA.")

    except Exception as e:

        print(f"\nERROR reading FASTA file for {ds}:")
        print(e)
        continue

    # =====================================================
    # READ COUNT TABLE
    # =====================================================

    try:

        df = pd.read_csv(
            count_file,
            sep="\t",
            skiprows=1
        )

        print(f"Loaded count table with {df.shape[0]} ASVs and {df.shape[1]-1} samples.")

    except Exception as e:

        print(f"\nERROR reading count table for {ds}:")
        print(e)
        continue

    # =====================================================
    # RENAME FIRST COLUMN
    # =====================================================

    df.rename(columns={df.columns[0]: "ASV_ID"}, inplace=True)

    # =====================================================
    # CLEAN SAMPLE NAMES
    # Remove '-nc'
    # =====================================================

    clean_columns = ["ASV_ID"]

    for col in df.columns[1:]:
        clean_columns.append(col.replace("-nc", ""))

    df.columns = clean_columns

    # =====================================================
    # REPLACE ASV IDs WITH REAL SEQUENCES
    # =====================================================

    original_rows = len(df)

    df["ASV_ID"] = df["ASV_ID"].map(seq_dict)

    # Remove rows where sequence not found
    df = df.dropna(subset=["ASV_ID"])

    removed_rows = original_rows - len(df)

    if removed_rows > 0:
        print(f"WARNING: Removed {removed_rows} ASVs with missing sequences.")

    # =====================================================
    # CONVERT COUNTS TO INTEGERS
    # =====================================================

    try:

        for col in df.columns[1:]:
            df[col] = df[col].astype(float).astype(int)

    except Exception as e:

        print(f"\nERROR converting counts to integers for {ds}:")
        print(e)
        continue

    # =====================================================
    # SET SEQUENCE COLUMN AS INDEX
    # =====================================================

    df = df.set_index("ASV_ID")

    # =====================================================
    # SAVE TIDIED COUNT TABLE
    # =====================================================

    try:

        df.to_csv(
            output_file,
            sep="\t"
        )

        print(f"\nTidied ASV table written to:")
        print(output_file)

    except Exception as e:

        print(f"\nERROR writing ASV table for {ds}:")
        print(e)
        continue

    # =====================================================
    # CREATE CLEAN FASTA FILE
    # =====================================================

    try:

        written_count = 0

        with open(clean_fasta_output, "w") as out_fasta:

            for original_id, sequence in seq_dict.items():

                clean_id = None

                # =================================================
                # FORMAT 1
                # Example:
                # SRR11949832-nc|SRR11949832.100190_ex
                # -> ASV_100190
                # =================================================

                match_numeric = re.search(r"\.(\d+)_ex$", original_id)

                if match_numeric:

                    asv_number = match_numeric.group(1)
                    clean_id = f"ASV_{asv_number}"

                else:

                    # =============================================
                    # FORMAT 2
                    # Example:
                    # WOF37-nc|1_1104_12962_4886_1_ex
                    # -> ASV_1_1104_12962_4886_1
                    # =============================================

                    match_complex = re.search(r"\|(.+)_ex$", original_id)

                    if match_complex:

                        asv_number = match_complex.group(1)
                        clean_id = f"ASV_{asv_number}"

                # =============================================
                # WRITE OUTPUT
                # =============================================

                if clean_id is not None:

                    out_fasta.write(f">{clean_id}\n")
                    out_fasta.write(f"{sequence}\n")

                    written_count += 1

                else:

                    print("WARNING: Could not parse FASTA header:")
                    print(original_id)

        print(f"\nClean FASTA written to:")
        print(clean_fasta_output)

        print(f"Sequences written: {written_count}")

    except Exception as e:

        print(f"\nERROR writing clean FASTA for {ds}:")
        print(e)
        continue

print("\nAll datasets processed.")