import pandas as pd
import os
import re
from Bio import SeqIO
import numpy as np

#===================================================================
# set variables -----------------------------------------------------
#===================================================================

datasets = ["Windermere_2017", "Marchamley", "RingTrial_Sean"] # folder names
base_path = ["Data/Processed"]
denoise_method = ["DADA2", "VSEARCH", "MetaBEAT"] # denoising methods
taxonomy_method = ["RDP", "MAPseq", "BLAST"]
base_path = "/Workspace/Users/dina.simons@environment-agency.gov.uk/Fish_pipeline_comparisons/Data/Processed/"

#====================================================================
# build ASV master table --------------------------------------------
#====================================================================

print("--------------------------------------")
print("\nImporting and merging ASVs...\n")
print("--------------------------------------")

master_list = []

for method in denoise_method:
    for ds in datasets:
 
        print(f"\nProcessing ASVs: {ds} | {method}")

        dataset_name = ds.lower()
        method_name = method.lower()

        count_file = os.path.join(base_path, ds, f"06_ASV_counts_{method}.tsv")
        fasta_file = os.path.join(base_path, ds, f"06_ASV_seqs_{method}.fasta")

        print(f"\nCount file location: {count_file}")
        print(f"Fasta file location: {fasta_file}")

        if not os.path.exists(count_file) or not os.path.exists(fasta_file):
            print(f"Missing files for {ds} | {method}")
            continue

        # -------------------------
        # READ COUNT TABLE
        # -------------------------
        asv_count = pd.read_csv(
            count_file,
            sep="\t",
            index_col=0
        )

        asv_count.index = asv_count.index.str.upper().str.strip()
        asv_count["sequence"] = asv_count.index

        # -------------------------
        # READ FASTA
        # -------------------------
        fasta_records = [
            (rec.id, str(rec.seq).upper().strip())
            for rec in SeqIO.parse(fasta_file, "fasta")
        ]

        asv_fasta_df = pd.DataFrame(fasta_records, columns=["asv", "sequence"])

        # -------------------------
        # RESHAPE (wide → long)
        # -------------------------
        count_long = (
            asv_count
            .reset_index(drop=True)
            .melt(id_vars="sequence", var_name="sampleID", value_name="reads")
        )

        count_long = count_long[count_long["reads"] > 0]

        # -------------------------
        # MERGE
        # -------------------------
        merged = pd.merge(count_long, asv_fasta_df, on="sequence", how="left")

        if merged.empty:
            print("Merge failed")
            continue

        merged["denoise_method"] = method_name
        merged["dataset"] = dataset_name

        merged = merged[["asv", "sampleID", "reads", "sequence", "denoise_method", "dataset"]]
        merged.columns = map(str.lower, merged.columns)

        master_list.append(merged)

# Combine
ASV_merged_dataset = pd.concat(master_list, ignore_index=True)
ASV_merged_dataset = ASV_merged_dataset.drop_duplicates()

# Convert to category (factors)
for col in ["asv", "sampleid", "denoise_method", "dataset"]:
    ASV_merged_dataset[col] = ASV_merged_dataset[col].astype("category")

print("\nRows of master ASV file:", len(ASV_merged_dataset))

#====================================================================
# taxonomy master table ---------------------------------------------
#====================================================================

print("--------------------------------------")
print("\nImporting and merging taxonomy...\n")
print("--------------------------------------")

def make_taxonomy_string(df):
    return df[["kingdom","phylum","class","order","family","genus","species"]]\
        .apply(lambda x: ";".join(x.fillna("unclassified").astype(str)), axis=1)
    
# label mapseq columns
mapseq_cols = [
    "asv", "hit_id", "bitscore", "identity", "matches", "mismatches", "gaps",
    "query_start", "query_end", "dbhit_start", "dbhit_end", "strand", "kingdom", "kingdom_combined_cf", "kingdom_score_cf",
    "phylum", "phylum_combined_cf", "phylum_score_cf",
    "class", "class_combined_cf", "class_score_cf",
    "order", "order_combined_cf", "order_score_cf",
    "family", "family_combined_cf", "family_score_cf",
    "genus", "genus_combined_cf", "genus_score_cf",
    "species", "species_combined_cf", "species_score_cf"
]

# loop

taxonomy_master_list = []
taxonomy_rdp_only_list = []

for ds in datasets:
    for method in denoise_method:

        dataset_name = ds.lower()
        method_name = method.lower()

        folder_path = os.path.join(base_path, ds)

        print("\nProcessing taxa:", dataset_name, method_name)

        files = os.listdir(folder_path)

        # MAPSEQ
        mapseq_file = [
            f for f in files if "mapseq" in f.lower()
            and "riaz" in f.lower()
            and f"_{method_name}_" in f.lower()
        ]

        if len(mapseq_file) == 1:
            print(f"\nMAPseq file location: {mapseq_file}")

            df_map = pd.read_csv(os.path.join(folder_path, mapseq_file[0]), sep="\t", header=0)
            print("Rows:", len(df_map))

            df_map = df_map.dropna(axis=1, how="all")
            df_map.columns = mapseq_cols
            df_map.columns = map(str.lower, df_map.columns)

            tax_map = df_map[["asv","kingdom","phylum","class","order","family","genus","species"]].copy()
            tax_map["taxonomy"] = make_taxonomy_string(tax_map)
            tax_map["taxonomy_method"] = "mapseq"
            tax_map["dataset"] = dataset_name
            tax_map["denoise_method"] = method_name

        taxonomy_master_list.append(tax_map)

        # RDP
        rdp_file = [f for f in files if "assigned_taxonomy" in f.lower() 
                    and "riaz" in f.lower()
                    and method_name in f.lower()]

        if len(rdp_file) == 1:
            print(f"\nRDP file location: {rdp_file}")

            df_rdp = pd.read_csv(os.path.join(folder_path, rdp_file[0]))
            print("Rows:", len(df_rdp))

            df_rdp = df_rdp.reset_index().rename(columns={"index": "sequence"})
            df_rdp.columns = map(str.lower, df_rdp.columns)

            tax_rdp = df_rdp[["sequence","kingdom","phylum","class","order","family","genus","species"]].copy()
            tax_rdp["taxonomy"] = make_taxonomy_string(tax_rdp)
            tax_rdp["taxonomy_method"] = "RDP"
            tax_rdp["dataset"] = dataset_name
            tax_rdp["denoise_method"] = method_name

        taxonomy_rdp_only_list.append(tax_rdp)

        # LCA
        lca_file = [f for f in files if "lca" in f.lower() 
                    and "riaz" in f.lower()
                    and method_name in f.lower()]

        if len(lca_file) == 1:
            print(f"\nBLAST LCA file location: {lca_file}")

            df_lca = pd.read_csv(os.path.join(folder_path, lca_file[0]), sep="\t")
            print("Rows:", len(df_lca))

            df_lca["kingdom"] = np.nan
            df_lca = df_lca.rename(columns={"ASV_ID": "asv"})
            df_lca.columns = map(str.lower, df_lca.columns)

            tax_lca = df_lca[["asv","kingdom","phylum","class","order","family","genus","species"]].copy()
            tax_lca["taxonomy"] = make_taxonomy_string(tax_lca)
            tax_lca["taxonomy_method"] = "BLAST"
            tax_lca["dataset"] = dataset_name
            tax_lca["denoise_method"] = method_name

        taxonomy_master_list.append(tax_lca)

#====================================================================
# combine and merge -------------------------------------------------
#====================================================================

print("--------------------------------------")
print("\nMerging and summarising ASVs and taxa...\n")
print("--------------------------------------")

taxonomy_master_dataset = pd.concat(taxonomy_master_list)
taxonomy_rdp_dataset = pd.concat(taxonomy_rdp_only_list) # because rdp only provides seqs

# remove duplicates
taxonomy_master_dataset = taxonomy_master_dataset.drop_duplicates(
    subset=["asv","taxonomy_method","dataset","denoise_method"]
)

taxonomy_rdp_dataset = taxonomy_rdp_dataset.drop_duplicates(
    subset=["sequence","taxonomy_method","dataset","denoise_method"]
)

# merge
# rdp first by sequence
asv_rdp_join = pd.merge(ASV_merged_dataset,
    taxonomy_rdp_dataset,
    on=["sequence","dataset","denoise_method"],
    how="left"
)

print("\nFinal RDP rows:", len(asv_rdp_join))
print("Missing taxonomy:", asv_rdp_join["taxonomy"].isna().sum())

asv_other_join = pd.merge(ASV_merged_dataset,
    taxonomy_master_dataset,
    on=["asv","dataset", "denoise_method"],
    how="left"
)

asv_other_join = asv_other_join[
    asv_other_join["taxonomy_method"].notna()
]

print("\nFinal taxa for mapseq and blast rows:", len(asv_other_join))
print("Missing taxonomy:", asv_other_join["taxonomy"].isna().sum())

# join the two
final_master_dataset = pd.concat([asv_rdp_join, asv_other_join], ignore_index=True)
final_master_dataset = final_master_dataset.drop_duplicates()

#final_master_dataset = final_master_dataset[
    #final_master_dataset["taxonomy_method"].notna()
#]

print("\nFinal total rows:", len(final_master_dataset))
print("Missing taxonomy:", final_master_dataset["taxonomy"].isna().sum())

# save
final_master_dataset.to_csv("Data/Processed/master_long_data.csv", index=False)