import pandas as pd
from pathlib import Path

MANUAL_DIR = Path(__file__).parent
DATASETS = ["GSE21138", "GSE38485", "GSE38484", "GSE48072"]
TOP_N = 5000

for dataset in DATASETS:
    print(f"[INFO] Processing {dataset}...")

    ds_dir = MANUAL_DIR / dataset
    de_file = ds_dir / "differential_expression.tsv"
    processed_file = ds_dir / "processed_expression.tsv"

    de_df = pd.read_csv(de_file, sep="\t")

    required_cols = ["gene_symbol", "logFC", "p_value"]
    if not all(col in de_df.columns for col in required_cols):
        raise ValueError(f"{dataset}: differential_expression.tsv missing required columns {required_cols}")

    de_df["ranking_score"] = de_df["logFC"].abs()

    top_genes = de_df.sort_values("ranking_score", ascending=False).head(TOP_N)

    top_genes_file = ds_dir / "top_genes.tsv"
    top_genes.to_csv(top_genes_file, sep="\t", index=False)
    print(f"[INFO] Saved {top_genes_file}")

    top_genes_mapped = top_genes.groupby("gene_symbol", as_index=False).agg({
        "logFC": "max",
        "p_value": "min",
        "ranking_score": "max"
    })
    top_genes_mapped_file = ds_dir / "top_genes_mapped.tsv"
    top_genes_mapped.to_csv(top_genes_mapped_file, sep="\t", index=False)
    print(f"[INFO] Saved {top_genes_mapped_file}")

print("[INFO] All datasets processed.")
