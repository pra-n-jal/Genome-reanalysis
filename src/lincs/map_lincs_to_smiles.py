
from pathlib import Path
import re
import argparse
import logging
import pandas as pd
import numpy as np

logging.basicConfig(level=logging.INFO, format='[%(levelname)s] %(message)s')

CANDIDATE_PERT_FILES = [
    "GSE92742_Broad_LINCS_pert_info.txt",
    "GSE92742_Broad_LINCS_pert_info.txt.gz",
    "pert_info.txt",
    "perturbagen_info.txt",
    "pert_info.tsv",
    "pert_info.csv",
    "perturbagen_metadata.tsv",
    "*pert_info*.txt",
    "*pert_info*.tsv",
    "*pert*info*.txt",
]

SMILES_COL_CANDIDATES = ["canonical_smiles", "canonical_smiles_smiles", "smiles", "SMILES", "chembl_smiles", "inchi_key"]
PERT_ID_COL_CANDIDATES = ["pert_id", "pert_id", "pert_iname", "perturbagen_id", "perturbagen", "pert_id"]

def find_pert_file(search_root: Path) -> Path:
    logging.info("Searching for perturbagen metadata file...")
    for candidate in CANDIDATE_PERT_FILES:
        if "*" in candidate:
            matched = list(search_root.glob(candidate))
            if matched:
                logging.info(f"Found pert file (glob): {matched[0].name}")
                return matched[0]
        else:
            p = search_root / candidate
            if p.exists():
                logging.info(f"Found pert file: {p.name}")
                return p

    for p in search_root.iterdir():
        if p.is_file() and ('pert' in p.name.lower() and 'info' in p.name.lower()):
            logging.info(f"Found fallback pert file: {p.name}")
            return p
    logging.warning("No perturbagen metadata file found automatically.")
    return None

def load_pert_info(path: Path) -> pd.DataFrame:
    if path is None:
        return None
    try:
        df = pd.read_csv(path, sep="\t", low_memory=False)
    except Exception:
        try:
            df = pd.read_csv(path, sep=",", low_memory=False)
        except Exception as e:
            logging.error(f"Could not read pert info file {path}: {e}")
            return None

    df.columns = [c.strip() for c in df.columns]
    return df

def detect_cols(df: pd.DataFrame):
    cols = list(df.columns)
    lowermap = {c.lower(): c for c in cols}
    smiles_col = None
    id_col = None
    name_col = None
    for cand in SMILES_COL_CANDIDATES:
        if cand in cols:
            smiles_col = cand
            break
        if cand.lower() in lowermap:
            smiles_col = lowermap[cand.lower()]
            break
    for cand in PERT_ID_COL_CANDIDATES:
        if cand in cols:
            id_col = cand
            break
        if cand.lower() in lowermap:
            id_col = lowermap[cand.lower()]
            break

    for candidate in ("pert_iname","pert_iname","pert_name","name","compound_name"):
        if candidate in cols:
            name_col = candidate
            break
        if candidate.lower() in lowermap:
            name_col = lowermap[candidate.lower()]
            break
    return id_col, smiles_col, name_col

def extract_brd_ids(series: pd.Series) -> pd.Series:

    brd = series.astype(str).str.extract(r'(BRD-[A-Za-z0-9\-\_]+)', expand=False)

    def fallback_extract(s):
        if not isinstance(s, str):
            return None

        toks = re.split(r'[:;|,\s]', s)
        for t in toks:
            if t.startswith('BRD-'):
                return t
        m = re.search(r'(BRD-[A-Za-z0-9\-\_]+)', s)
        return m.group(1) if m else None
    brd = brd.where(brd.notna(), series.apply(fallback_extract))

    brd = brd.astype(str).str.strip().replace({'nan': None})
    brd = brd.where(brd != 'None', None)
    return brd

def main(args):
    root = Path('.')
    labeled_path = Path(args.labeled)
    if not labeled_path.exists():
        logging.error(f"Labeled file not found: {labeled_path}")
        return

    logging.info(f"Loading labeled file: {labeled_path}")
    df = pd.read_csv(labeled_path, sep="\t", low_memory=False)
    if 'compound' not in df.columns:
        logging.error("Expected 'compound' column in labeled file.")
        return

    df['compound'] = df['compound'].astype(str)
    df['BRD_ID'] = extract_brd_ids(df['compound'])
    n_brd = df['BRD_ID'].notna().sum()
    logging.info(f"Extracted BRD IDs for {n_brd} / {len(df)} rows")

    pert_file = find_pert_file(root)
    pert_df = load_pert_info(pert_file) if pert_file else None
    if pert_df is None:

        possibles = list(Path.cwd().glob('**/*pert_info*.txt')) + list(Path.home().glob('**/*pert_info*.txt'))
        if possibles:
            pert_df = load_pert_info(possibles[0])
            if pert_df is not None:
                logging.info(f"Loaded pert_info from fallback: {possibles[0].name}")
    if pert_df is None:
        logging.error("No perturbagen metadata found. Place 'GSE92742_Broad_LINCS_pert_info.txt' or similar in working dir.")
        return

    id_col, smiles_col, name_col = detect_cols(pert_df)
    logging.info(f"Detected pert_info columns -> id: {id_col}, smiles: {smiles_col}, name: {name_col}")

    if id_col:
        pert_df[id_col] = pert_df[id_col].astype(str).str.strip()

    merged = df.copy()
    if id_col and smiles_col and id_col in pert_df.columns and smiles_col in pert_df.columns:

        use_cols = [id_col, smiles_col]
        if name_col and name_col in pert_df.columns:
            use_cols.append(name_col)
        pert_sub = pert_df[use_cols].copy()

        merged = merged.merge(pert_sub, left_on='BRD_ID', right_on=id_col, how='left')

        if smiles_col in merged.columns:
            merged.rename(columns={smiles_col: 'smiles'}, inplace=True)
        if name_col and name_col in merged.columns:
            merged.rename(columns={name_col: 'pert_iname'}, inplace=True)
    else:

        if 'compound' in df.columns and name_col and name_col in pert_df.columns:
            logging.info("Falling back to merge on pert_iname / compound name")
            pert_sub = pert_df[[name_col]].copy()
            pert_sub = pert_sub.rename(columns={name_col: 'pert_iname'})
            merged = merged.merge(pert_sub, left_on='compound', right_on='pert_iname', how='left')
        else:
            logging.warning("Could not perform BRD->SMILES merge: required columns missing in pert_info.")

            merged = merged

    has_smiles = 'smiles' in merged.columns and merged['smiles'].notna().sum() > 0
    n_mapped = merged['smiles'].notna().sum() if 'smiles' in merged.columns else 0
    logging.info(f"Mapped {n_mapped} / {len(merged)} rows to SMILES")

    if 'spearman_rank_pct' in merged.columns:

        merged['score'] = -pd.to_numeric(merged['spearman_rank_pct'], errors='coerce')
    elif 'spearman' in merged.columns:
        merged['score'] = -pd.to_numeric(merged['spearman'], errors='coerce')
    elif 'cosine' in merged.columns:
        merged['score'] = pd.to_numeric(merged['cosine'], errors='coerce')
    else:
        merged['score'] = 0.0

    pct = float(args.percentile)
    if merged['score'].notna().sum() == 0:
        logging.warning("No numeric scores found; label_binary will be 0 for all rows.")
        merged['label_binary'] = 0
    else:
        cutoff = merged['score'].quantile(1.0 - (pct/100.0))
        merged['label_binary'] = (merged['score'] >= cutoff).astype(int)

    final_cols = ['smiles','score','label_binary','spearman','cosine','spearman_rank_pct','BRD_ID','pert_iname']
    final_present = [c for c in final_cols if c in merged.columns]
    final_df = merged[final_present].copy()

    if 'smiles' not in final_df.columns and 'pert_iname' in final_df.columns:
        final_df = final_df.rename(columns={'pert_iname': 'compound_name'})

    outpath = Path("C4A_training_dataset.tsv")
    final_df.to_csv(outpath, sep="\t", index=False)
    logging.info(f"Wrote ML dataset to {outpath} ({len(final_df)} rows)")

    unmapped = merged.loc[merged['smiles'].isna() if 'smiles' in merged.columns else merged['BRD_ID'].isna(), 'BRD_ID']
    unmapped = unmapped.dropna().unique().tolist()
    unmapped_path = Path("unmapped_brd_ids.txt")
    unmapped_path.write_text("\n".join([str(x) for x in unmapped]))
    logging.info(f"Wrote {len(unmapped)} unmapped BRD IDs to {unmapped_path}")

if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--labeled', default='labeled_compound_C4A.tsv', help='Labeled compounds TSV from previous step')
    ap.add_argument('--percentile', type=float, default=5.0, help='Top percentile to call label_binary=1 (default 5% = top 5%)')
    args = ap.parse_args()
    main(args)
