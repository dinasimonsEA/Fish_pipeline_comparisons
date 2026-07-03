# script to tidy and filter mapseq output
# new version of mapseq doesn't allow for the simple format, so we need to filter manually by combined score columns
# format similar to other taxonomy assignment outputs

from pathlib import Path
import pandas as pd

THRESHOLD = 0.5

def process_mapseq(file_path):

    df = pd.read_csv(file_path, sep="\t")

    cf_cols = [c for c in df.columns if str(c).startswith("combined_cf")]

    rename_dict = {
        cf_cols[0]: "kingdom_cf",
        cf_cols[1]: "phylum_cf",
        cf_cols[2]: "class_cf",
        cf_cols[3]: "order_cf",
        cf_cols[4]: "family_cf",
        cf_cols[5]: "genus_cf",
        cf_cols[6]: "species_cf",
    }

    df = df.rename(columns=rename_dict)

    out = pd.DataFrame({
        "ASV_ID": df["#query"],
        "kingdom": df["MetaFish:Kingdom"].where(df["kingdom_cf"] >= THRESHOLD, "NA"),
        "phylum": df["Phylum"].where(df["phylum_cf"] >= THRESHOLD, "NA"),
        "class": df["Class"].where(df["class_cf"] >= THRESHOLD, "NA"),
        "order": df["Order"].where(df["order_cf"] >= THRESHOLD, "NA"),
        "family": df["Family"].where(df["family_cf"] >= THRESHOLD, "NA"),
        "genus": df["Genus"].where(df["genus_cf"] >= THRESHOLD, "NA"),
        "species": df["Species"].where(df["species_cf"] >= THRESHOLD, "NA")
    })

    return out


files = list(Path("Data/Processed").rglob("*.mseq"))

for f in files:

    result = process_mapseq(f)

    out_file = f.with_name(
        f.stem + "_taxonomy_cf05.tsv"
    )

    result.to_csv(
        out_file,
        sep="\t",
        index=False,
        na_rep="NA"
    )

    print(f"Saved {out_file}")