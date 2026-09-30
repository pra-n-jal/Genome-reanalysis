
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import re
import sys

TARGET_GENES = ["C4A", "TCF4", "CACNA1C", "CHD8"]

BASE = Path(".")
OUT_SUMMARY = BASE / "core_gene_summary.tsv"
OUT_PRESENCE = BASE / "core_gene_presence.tsv"
OUT_HEATMAP = BASE / "core_gene_presence_heatmap.png"

CANDIDATE_FILES = ["top_genes_mapped.tsv",
                   "top_genes_mapped_annotated_collapsed.tsv",
                   "top_genes_mapped_annotated.tsv",
                   "top_genes.tsv",
                   "processed_expression_top_genes.tsv"]

def find_gene_column(df):

    candidates = [c for c in df.columns if re.search(r'(^gene$|gene_symbol|symbol|hgnc|entrez)', c, re.I)]
    if candidates:
        return candidates[0]

    return None

def find_col_like(df, patterns):
    for p in patterns:
        for c in df.columns:
            if re.search(p, c, re.I):
                return c
    return None

summary_rows = []
presence = {}

for gse in sorted(BASE.glob("GSE*")):

    found = None
    for fname in CANDIDATE_FILES:
        f = gse / fname
        if f.exists():
            found = f
            break
    if not found:
        print(f"[WARN] {gse.name}: no candidate top_genes file found, skipping")
        continue

    try:
        df = pd.read_csv(found, sep="\t", low_memory=False)
    except Exception as e:
        print(f"[ERROR] Failed reading {found}: {e}", file=sys.stderr)
        continue

    gene_col = find_gene_column(df)
    if gene_col is None:

        gene_col = df.columns[0]

        if pd.api.types.is_numeric_dtype(df[gene_col]):
            print(f"[WARN] {gse.name}: guessed gene column {gene_col} looks numeric. File: {found}")

    logfc_col = find_col_like(df, [r'log2?fc', r'log_fc', r'logfc', r'log2.fold', r'foldchange', r'fold_change'])
    ranking_col = find_col_like(df, [r'ranking_score', r'rank(?!\w)|ranking', r'score'])
    adjp_col = find_col_like(df, [r'adj_p', r'padj', r'fdr', r'adj\.p', r'adjP', r'adjusted'])

    if gene_col in df.columns:
        df['__gene_up'] = df[gene_col].astype(str).str.upper().str.replace(r'\s+', '', regex=True)
    else:
        df['__gene_up'] = df.index.astype(str).str.upper().str.replace(r'\s+', '', regex=True)

    presence[gse.name] = {}
    for gene in TARGET_GENES:
        g_up = gene.upper().replace(" ", "")
        mask = df['__gene_up'] == g_up
        if mask.any():
            row = df.loc[mask].iloc[0]
            logfc = row[logfc_col] if logfc_col in row.index else ""
            ranking = row[ranking_col] if ranking_col in row.index else ""
            adjp = row[adjp_col] if adjp_col in row.index else ""
            summary_rows.append({
                "dataset": gse.name,
                "gene": gene,
                "present": 1,
                "logFC": logfc,
                "ranking_score": ranking,
                "adj_p": adjp,
                "source_file": found.name
            })
            presence[gse.name][gene] = 1
        else:
            summary_rows.append({
                "dataset": gse.name,
                "gene": gene,
                "present": 0,
                "logFC": "",
                "ranking_score": "",
                "adj_p": "",
                "source_file": found.name
            })
            presence[gse.name][gene] = 0

summary_df = pd.DataFrame(summary_rows)
summary_df.to_csv(OUT_SUMMARY, sep="\t", index=False)
print(f"[OK] wrote {OUT_SUMMARY}")

if presence:
    pres_df = pd.DataFrame.from_dict(presence, orient='index').fillna(0).astype(int)

    pres_df = pres_df[TARGET_GENES]
    pres_df = pres_df.transpose()
    pres_df.to_csv(OUT_PRESENCE, sep="\t")
    print(f"[OK] wrote {OUT_PRESENCE}")

    plt.figure(figsize=(max(4, len(pres_df.columns)*0.6), max(2, len(pres_df.index)*0.6)))

    plt.imshow(pres_df.values, aspect='auto', interpolation='nearest')
    plt.yticks(range(len(pres_df.index)), pres_df.index)
    plt.xticks(range(len(pres_df.columns)), pres_df.columns, rotation=45, ha='right')
    plt.title("Presence (1) / Absence (0) of target genes in top_genes_mapped")
    plt.colorbar(label='presence')
    plt.tight_layout()
    plt.savefig(OUT_HEATMAP, dpi=200)
    plt.close()
    print(f"[OK] wrote {OUT_HEATMAP}")
else:
    print("[WARN] No GSE folders processed. Check directory or filenames.")
