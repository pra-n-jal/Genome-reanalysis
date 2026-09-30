import os
import pandas as pd

parent_dir = r"C:\Users\pranj\scz_connectivity_multi\Final_results\Manual"

c4_genes = ["C4A", "C4B"]

rows = []

for gse_folder in os.listdir(parent_dir):
    folder_path = os.path.join(parent_dir, gse_folder)
    if os.path.isdir(folder_path):
        de_file = os.path.join(folder_path, "differential_expression.tsv")
        if os.path.exists(de_file):
            df = pd.read_csv(de_file, sep="\t")

            mask = df['gene_symbol'].str.contains('|'.join(c4_genes))
            c4_df = df[mask].copy()
            if not c4_df.empty:
                for _, row in c4_df.iterrows():
                    rows.append({
                        "dataset": gse_folder,
                        "gene_symbol": row['gene_symbol'],
                        "tstat": row.get('tstat', ''),
                        "p_value": row.get('p_value', ''),
                        "logFC": row.get('logFC', ''),
                        "adj_p": row.get('adj_p', '')
                    })

out_df = pd.DataFrame(rows)
out_file = os.path.join(parent_dir, "C4A_genes_combined.tsv")
out_df.to_csv(out_file, sep="\t", index=False)

print(f"[INFO] Saved C4A/B gene rows to {out_file}")
