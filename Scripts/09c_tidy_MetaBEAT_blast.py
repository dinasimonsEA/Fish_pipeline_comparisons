import pandas as pd
import ast
import os

# ====================================================
# DATASETS
# ====================================================

datasets = [
    "Windermere_2017",
    "RingTrial_Sean",
    "Marchamley"
]

# ====================================================
# LOOP THROUGH DATASETS
# ====================================================

for ds in datasets:

    print("\n" + "=" * 70)
    print(f"Processing: {ds}")
    print("=" * 70)

    base_path = f"Data/Processed/{ds}"

    taxonomy_file = (
        f"{base_path}/08_MetaBEAT_blast_95.tsv"
    )

    lookup_file = (
        f"{base_path}/06_OTU_lookup_MetaBEAT.tsv"
    )

    output_file = (
        f"{base_path}/09_MetaBEAT_OTU_taxonomy_blast_95.tsv"
    )

    # ====================================================
    # CHECK FILES EXIST
    # ====================================================

    if not os.path.exists(taxonomy_file):

        print(f"Missing taxonomy file:\n{taxonomy_file}")
        continue

    if not os.path.exists(lookup_file):

        print(f"Missing lookup file:\n{lookup_file}")
        continue

    try:

        # ====================================================
        # READ LOOKUP TABLE
        # ====================================================

        lookup = pd.read_csv(
            lookup_file,
            sep="\t"
        )

        id_map = dict(
            zip(
                lookup["Original_OTU"],
                lookup["New_OTU"]
            )
        )

        # ====================================================
        # READ TAXONOMY FILE
        # ====================================================

        df = pd.read_csv(
            taxonomy_file,
            sep="\t",
            skiprows=1
        )

        df.rename(
            columns={
                df.columns[0]: "Original_OTU",
                df.columns[-1]: "taxonomy"
            },
            inplace=True
        )

        df = df[
            ["Original_OTU", "taxonomy"]
        ].copy()

        # ====================================================
        # APPLY NEW OTU IDS
        # ====================================================

        df["ASV_ID"] = df["Original_OTU"].map(id_map)

        missing_ids = df["ASV_ID"].isna().sum()

        if missing_ids > 0:

            print(
                f"WARNING: {missing_ids} IDs "
                f"could not be mapped."
            )

        # ====================================================
        # TAXONOMY PARSER
        # ====================================================

        def parse_taxonomy(tax_string):

            ranks = {
                "Species": "NA",
                "Genus": "NA",
                "Family": "NA",
                "Order": "NA",
                "Class": "NA",
                "Phylum": "NA",
                "Kingdom": "NA",
                "lca_rank": "NA",
                "lca_name": "NA"
            }

            if pd.isna(tax_string):
                return ranks

            try:

                taxa = ast.literal_eval(tax_string)

                for item in taxa:

                    if item.startswith("k__"):
                        ranks["Kingdom"] = item[3:]

                    elif item.startswith("p__"):
                        ranks["Phylum"] = item[3:]

                    elif item.startswith("c__"):
                        ranks["Class"] = item[3:]

                    elif item.startswith("o__"):
                        ranks["Order"] = item[3:]

                    elif item.startswith("f__"):
                        ranks["Family"] = item[3:]

                    elif item.startswith("g__"):
                        ranks["Genus"] = item[3:]

                    elif item.startswith("s__"):
                        ranks["Species"] = (
                            item[3:]
                            .replace("_", " ")
                        )

                hierarchy = [
                    ("Species", ranks["Species"]),
                    ("Genus", ranks["Genus"]),
                    ("Family", ranks["Family"]),
                    ("Order", ranks["Order"]),
                    ("Class", ranks["Class"]),
                    ("Phylum", ranks["Phylum"]),
                    ("Kingdom", ranks["Kingdom"])
                ]

                for rank, value in hierarchy:

                    if value != "NA":

                        ranks["lca_rank"] = rank
                        ranks["lca_name"] = value

                        break

            except Exception:
                pass

            return ranks

        # ====================================================
        # PARSE TAXONOMY
        # ====================================================

        tax_df = pd.DataFrame(
            df["taxonomy"]
            .apply(parse_taxonomy)
            .tolist()
        )

        # ====================================================
        # CREATE OUTPUT TABLE
        # ====================================================

        final_df = pd.concat(
            [
                df["ASV_ID"],
                tax_df
            ],
            axis=1
        )

        final_df = final_df[
            [
                "ASV_ID",
                "Species",
                "Genus",
                "Family",
                "Order",
                "Class",
                "Phylum",
                "Kingdom",
                "lca_rank",
                "lca_name"
            ]
        ]

        final_df = (
            final_df
            .dropna(subset=["ASV_ID"])
            .drop_duplicates(subset="ASV_ID")
            .sort_values("ASV_ID")
        )

        # ====================================================
        # WRITE OUTPUT
        # ====================================================

        final_df.to_csv(
            output_file,
            sep="\t",
            index=False
        )

        print(
            f"Written: {output_file}"
        )
        print(
            f"Rows: {len(final_df):,}"
        )

    except Exception as e:

        print(
            f"ERROR processing {ds}:"
        )
        print(e)

        continue

print("\nAll available datasets processed.")