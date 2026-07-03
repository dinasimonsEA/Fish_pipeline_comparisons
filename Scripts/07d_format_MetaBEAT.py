import pandas as pd
from collections import OrderedDict
import re
import os

# =========================================================
# DATASETS
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

    print("\n" + "=" * 70)
    print(f"Processing dataset: {ds}")
    print("=" * 70)

    # =====================================================
    # FILE PATHS
    # =====================================================

    base_path = (
        f"/Workspace/Users/dina.simons@environment-agency.gov.uk/"
        f"Fish_pipeline_comparisons/Data/Processed/{ds}"
    )

    count_file = f"{base_path}/06_ASV_counts_MetaBEAT_untidied.tsv"

    fasta_file = (
        f"{base_path}/06_ASV_seqs_MetaBEAT_untidied.fasta"
    )

    output_count = (
        f"{base_path}/06_ASV_counts_MetaBEAT.tsv"
    )

    output_fasta = (
        f"{base_path}/06_ASV_seqs_MetaBEAT.fasta"
    )

    output_lookup = (
    f"{base_path}/06_OTU_lookup_MetaBEAT.tsv"
    )

    # =====================================================
    # CHECK FILES EXIST
    # =====================================================

    if not os.path.exists(count_file):

        print(f"\nMissing count file:\n{count_file}")
        continue

    if not os.path.exists(fasta_file):

        print(f"\nMissing FASTA file:\n{fasta_file}")
        continue

    # =====================================================
    # READ FASTA
    # =====================================================

    raw_to_sequence = OrderedDict()
    otu_to_sequence = OrderedDict()
    sequence_to_otu = OrderedDict()
    raw_to_new_otu = OrderedDict()

    otu_counter = 1

    try:

        with open(fasta_file, "r") as f:

            current_header = None
            current_seq = []

            for line in f:

                line = line.strip()

                # =============================================
                # FASTA HEADER
                # =============================================

                if line.startswith(">"):

                    # Save previous sequence
                    if current_header is not None:

                        sequence = "".join(current_seq)

                        # Reuse OTU if sequence already exists
                        if sequence in sequence_to_otu:

                            otu_id = sequence_to_otu[sequence]

                        else:

                            otu_id = f"OTU_{otu_counter:05d}"

                            sequence_to_otu[sequence] = otu_id
                            otu_to_sequence[otu_id] = sequence

                            otu_counter += 1

                        raw_to_new_otu[current_header] = otu_id
                        raw_to_sequence[current_header] = sequence

                    # Start new record
                    current_header = line[1:]
                    current_seq = []

                else:

                    current_seq.append(line)

            # =============================================
            # SAVE FINAL RECORD
            # =============================================

            if current_header is not None:

                sequence = "".join(current_seq)

                if sequence in sequence_to_otu:

                    otu_id = sequence_to_otu[sequence]

                else:

                    otu_id = f"OTU_{otu_counter:05d}"

                    sequence_to_otu[sequence] = otu_id
                    otu_to_sequence[otu_id] = sequence

                raw_to_new_otu[current_header] = otu_id
                raw_to_sequence[current_header] = sequence

            # save look-up table
            try:

                lookup_df = pd.DataFrame({
                    "Original_OTU": list(raw_to_new_otu.keys()),
                    "New_OTU": list(raw_to_new_otu.values())
                })

                lookup_df.to_csv(
                    output_lookup,
                    sep="\t",
                    index=False
                )

                print("\nLookup table written:")
                print(output_lookup)

            except Exception as e:

                print("\nERROR writing lookup table:")
                print(e)

        print(f"\nUnique OTUs: {len(otu_to_sequence)}")

    except Exception as e:

        print("\nERROR reading FASTA:")
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

        print(
            f"\nLoaded count table:"
            f"\nRows    : {df.shape[0]}"
            f"\nSamples : {df.shape[1] - 1}"
        )

    except Exception as e:

        print("\nERROR reading count table:")
        print(e)

        continue

    # =====================================================
    # RENAME FIRST COLUMN
    # =====================================================

    df.rename(
        columns={df.columns[0]: "RAW_ID"},
        inplace=True
    )

    # =====================================================
    # CLEAN SAMPLE NAMES
    # =====================================================

    clean_cols = ["RAW_ID"]

    for col in df.columns[1:]:

        clean_col = re.sub(r"-nc$", "", col)

        clean_cols.append(clean_col)

    df.columns = clean_cols

    # =====================================================
    # MAP RAW IDs TO DNA SEQUENCES
    # =====================================================

    df["Sequence"] = df["RAW_ID"].map(raw_to_sequence)

    missing = df["Sequence"].isna().sum()

    if missing > 0:

        print(
            f"\nWARNING: {missing} rows "
            f"could not be mapped to sequences."
        )

    # Remove failed mappings
    df = df.dropna(subset=["Sequence"])

    # =====================================================
    # REMOVE RAW IDS
    # =====================================================

    df.drop(columns=["RAW_ID"], inplace=True)

    # =====================================================
    # MOVE SEQUENCE TO FIRST COLUMN
    # =====================================================

    cols = ["Sequence"] + [
        c for c in df.columns if c != "Sequence"
    ]

    df = df[cols]

    # =====================================================
    # CONVERT COUNTS TO INTEGERS
    # =====================================================

    for col in df.columns[1:]:

        df[col] = (
            pd.to_numeric(df[col], errors="coerce")
            .fillna(0)
            .astype(int)
        )

    # =====================================================
    # SET SEQUENCE AS INDEX
    # =====================================================

    df = df.set_index("Sequence")

    # =====================================================
    # SAVE COUNT TABLE
    # =====================================================

    try:

        df.to_csv(
            output_count,
            sep="\t"
        )

        print("\nCount table written:")
        print(output_count)

    except Exception as e:

        print("\nERROR writing count table:")
        print(e)

        continue

    # =====================================================
    # WRITE CLEAN FASTA
    # =====================================================

    try:

        with open(output_fasta, "w") as out_fasta:

            for otu_id, sequence in otu_to_sequence.items():

                out_fasta.write(f">{otu_id}\n")
                out_fasta.write(f"{sequence}\n")

        print("\nFASTA written:")
        print(output_fasta)

    except Exception as e:

        print("\nERROR writing FASTA:")
        print(e)

        continue

print("\nAll datasets processed successfully.")