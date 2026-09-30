import os
import pandas as pd
import requests
import json
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

BASE_DIR = r"C:\Users\pranj\scz_connectivity_multi\Final_results\Manual"
OUTPUT_FILE = "manual_drug_scores.tsv"
PVAL_THRESH = 0.05
TSTAT_THRESH = 2
GENE_CAP = 1000

MAX_WORKERS = 6
MAX_RETRIES = 3

ENRICHR_ADD_URL = "https://maayanlab.cloud/Enrichr/addList"
ENRICHR_ENRICH_URL = "https://maayanlab.cloud/Enrichr/enrich"

LIBRARIES = [
    "LINCS_L1000_Chem_Pert_up",
    "LINCS_L1000_Chem_Pert_down",
    "Drug_Perturbations_from_GEO_up",
    "Drug_Perturbations_from_GEO_down"
]

def extract_gene_sets(diff_expr_file):
    df = pd.read_csv(diff_expr_file, sep="\t")

    if "gene_symbol" in df.columns:
        gene_col = "gene_symbol"
    elif "gene" in df.columns:
        gene_col = "gene"
    else:
        raise ValueError(f"{diff_expr_file} missing 'gene_symbol' or 'gene' column")

    colmap = {c.lower(): c for c in df.columns}
    if "tstat" in colmap and "pval" in colmap:
        tstat_col, pval_col = colmap["tstat"], colmap["pval"]
    elif "combined_z" in colmap and "combined_p" in colmap:
        tstat_col, pval_col = colmap["combined_z"], colmap["combined_p"]
    elif "t" in colmap and "p" in colmap:
        tstat_col, pval_col = colmap["t"], colmap["p"]
    else:
        raise ValueError(f"{diff_expr_file} has unsupported columns")

    up = (
        df[(df[pval_col] < PVAL_THRESH) & (df[tstat_col] > TSTAT_THRESH)]
        .sort_values(by=tstat_col, ascending=False)
        .head(GENE_CAP)[gene_col]
        .dropna()
        .unique()
        .tolist()
    )
    down = (
        df[(df[pval_col] < PVAL_THRESH) & (df[tstat_col] < -TSTAT_THRESH)]
        .sort_values(by=tstat_col, ascending=True)
        .head(GENE_CAP)[gene_col]
        .dropna()
        .unique()
        .tolist()
    )

    print(f"  Using {len(up)} up, {len(down)} down genes (capped at {GENE_CAP})")
    return up, down

def enrichr_submit(gene_list, description=""):
    genes_str = "\n".join(gene_list)
    payload = {"list": (None, genes_str), "description": (None, description)}

    for attempt in range(MAX_RETRIES):
        try:
            response = requests.post(ENRICHR_ADD_URL, files=payload, timeout=30)
            if response.ok:
                return json.loads(response.text)["userListId"]
        except Exception:
            time.sleep(2 * (attempt + 1))
    raise Exception("Failed to submit gene list to Enrichr after retries")

def enrichr_get_results(user_list_id, lib):
    for attempt in range(MAX_RETRIES):
        try:
            url = f"{ENRICHR_ENRICH_URL}?userListId={user_list_id}&backgroundType={lib}"
            response = requests.get(url, timeout=30)
            if response.ok:
                return lib, response.json()
        except Exception:
            time.sleep(2 * (attempt + 1))
    print(f"Warning: Failed to fetch results for {lib}")
    return lib, None

def run_manual_pipeline(dataset_dir):
    diff_expr_file = os.path.join(dataset_dir, "differential_expression.tsv")
    if not os.path.exists(diff_expr_file):
        print(f"No differential_expression.tsv found in {dataset_dir}")
        return

    print(f"\n=== Processing {diff_expr_file} ===")
    up, down = extract_gene_sets(diff_expr_file)

    all_results = []

    for gene_list, tag in [(up, "up"), (down, "down")]:
        if not gene_list:
            continue

        try:
            user_list_id = enrichr_submit(gene_list, description=f"{os.path.basename(dataset_dir)}_{tag}")
        except Exception as e:
            print(f"  !!! Failed submitting {tag} genes: {e}")
            continue

        with ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
            futures = [executor.submit(enrichr_get_results, user_list_id, lib) for lib in LIBRARIES]
            for fut in as_completed(futures):
                lib, data = fut.result()
                if not data or lib not in data:
                    continue
                for row in data[lib]:
                    rank, drug, pval, zscore, combined_score = row[0], row[1], row[2], row[3], row[4]
                    all_results.append([drug, lib, tag, pval, zscore, combined_score, rank])

    if all_results:
        df_out = pd.DataFrame(
            all_results,
            columns=["drug", "library", "direction", "pval", "zscore", "combined_score", "rank"],
        )
        out_path = os.path.join(dataset_dir, OUTPUT_FILE)
        df_out.to_csv(out_path, sep="\t", index=False)
        print(f"  Saved drug lookup -> {out_path}")
    else:
        print(f"  !!! No enrichment results for {dataset_dir}")

def main():
    for gse in os.listdir(BASE_DIR):
        dataset_dir = os.path.join(BASE_DIR, gse)
        if os.path.isdir(dataset_dir):
            run_manual_pipeline(dataset_dir)

if __name__ == "__main__":
    main()
