
from pathlib import Path
import argparse
import logging
import pandas as pd
import numpy as np
from scipy.spatial.distance import cosine
from scipy.stats import spearmanr

logging.basicConfig(level=logging.INFO, format='[%(levelname)s] %(message)s')

def read_signature(path: Path):
    if not path.exists():
        logging.error(f"Signature file not found: {path}")
        return None
    df = pd.read_csv(path, sep='\t', index_col=0)

    if 'combined_z' in df.columns:
        sig = df['combined_z'].copy()
    elif 'mean_logFC' in df.columns:
        sig = df['mean_logFC'].copy()
    else:

        numeric_cols = [c for c in df.columns if pd.api.types.is_numeric_dtype(df[c])]
        if not numeric_cols:
            logging.error("No numeric column found in signature file.")
            return None
        sig = df[numeric_cols[0]].copy()
    sig.index = sig.index.astype(str).str.upper().str.strip()
    sig = sig.dropna()
    logging.info(f"Loaded signature with {len(sig)} genes (using '{sig.name}')")
    return sig

def find_manual_scores(root: Path):

    candidates = list(root.glob('**/manual_drug_scores.tsv'))
    if not candidates:
        return None

    logging.info(f"Found manual drug score file: {candidates[0]}")
    return candidates[0]

def read_lincs(path: Path):

    if not path.exists():
        logging.error(f"LINCS file not found: {path}")
        return None

    if path.suffix.lower() == '.gctx':
        logging.info(f"[GCTX] Parsing LINCS GCTX file: {path}")
        try:
            from cmapPy.pandasGEXpress.parse_gctx import parse
        except ImportError:
            logging.error("Missing cmapPy. Install it with: pip install cmapPy")
            return None

        gct = parse(str(path), convert_neg_666=True)
        df = gct.data_df
        df.index = df.index.astype(str).str.upper().str.strip()
        logging.info(f"[GCTX] Loaded matrix with shape {df.shape} (genes x profiles)")
        return df

    try:
        df = pd.read_csv(path, sep='\t', index_col=0)
        if df.shape[0] < df.shape[1]:
            df.index = df.index.astype(str).str.upper().str.strip()
            logging.info(f"[TSV] Loaded wide LINCS matrix: {df.shape}")
            return df
    except Exception as e:
        logging.debug(f"Wide TSV read failed: {e}")

    try:
        df2 = pd.read_csv(path, sep='\t')
        cols = [c.lower() for c in df2.columns]
        gene_col = None
        compound_col = None
        value_col = None
        for c in df2.columns:
            cl = c.lower()
            if cl in ('gene','genes','symbol','gene_symbol'):
                gene_col = c
            if cl in ('compound','pert_iname','perturbagen','sig_id'):
                compound_col = c
            if cl in ('value','zscore','logfc','score','signature'):
                value_col = c
        if gene_col and compound_col and value_col:
            wide = df2.pivot_table(index=gene_col, columns=compound_col, values=value_col, aggfunc='mean')
            wide.index = wide.index.astype(str).str.upper().str.strip()
            logging.info(f"[TSV-long] Pivoted long LINCS to wide: {wide.shape}")
            return wide
    except Exception as e:
        logging.debug(f"Long TSV read failed: {e}")

    logging.error("Could not parse LINCS file format.")
    return None

def compute_scores(sig: pd.Series, lincs_wide: pd.DataFrame, method='spearman'):

    genes = sig.index.intersection(lincs_wide.index)
    if len(genes) == 0:
        logging.error("No overlapping genes between signature and LINCS data.")
        return None
    logging.info(f"Using {len(genes)} overlapping genes for scoring")
    sig_vec = sig.loc[genes].astype(float).values
    results = []
    for comp in lincs_wide.columns:
        comp_vec = lincs_wide.loc[genes, comp].astype(float).values

        mask = ~np.isnan(comp_vec) & ~np.isnan(sig_vec)
        if mask.sum() < max(3, int(len(genes)*0.1)):

            results.append((comp, np.nan, np.nan))
            continue
        v1 = sig_vec[mask]
        v2 = comp_vec[mask]

        try:
            rho, _ = spearmanr(v1, v2)
        except Exception:
            rho = np.nan

        try:
            cos = 1 - cosine(v1, v2)
        except Exception:
            cos = np.nan
        results.append((comp, float(rho) if np.isfinite(rho) else np.nan, float(cos) if np.isfinite(cos) else np.nan))
    resdf = pd.DataFrame(results, columns=['compound','spearman','cosine']).set_index('compound')
    return resdf

def label_by_percentile(resdf: pd.DataFrame, percentile_cut=5.0):

    resdf = resdf.copy()
    resdf['spearman_rank_pct'] = resdf['spearman'].rank(pct=True, method='average') * 100.0

    resdf['label'] = 'nonreverser'
    mask = (resdf['spearman'].notna()) & (resdf['spearman_rank_pct'] <= percentile_cut)
    resdf.loc[mask, 'label'] = 'reverser'
    return resdf

def main(args):
    sig = read_signature(Path(args.signature))
    if sig is None:
        return

    topn = args.topn
    up = sig[sig > 0].sort_values(ascending=False).head(topn)
    down = sig[sig < 0].sort_values(ascending=True).head(topn)
    up.to_csv('up_signature.tsv', sep='\t', header=True)
    down.to_csv('down_signature.tsv', sep='\t', header=True)
    logging.info(f"Wrote up_signature.tsv ({len(up)}) and down_signature.tsv ({len(down)})")

    lincs_path = Path(args.lincs) if args.lincs else None
    lincs_df = None
    if lincs_path and lincs_path.exists():
        lincs_df = read_lincs(lincs_path)
    else:

        f = find_manual_scores(Path('.'))
        if f:
            try:
                md = pd.read_csv(f, sep='\t')

                cols = [c.lower() for c in md.columns]
                if 'gene' in cols or 'pert_iname' in cols or 'symbol' in cols:

                    gene_col = None
                    for c in md.columns:
                        if c.lower() in ('gene','symbol','genes'):
                            gene_col = c
                            break
                    comp_col = None
                    for c in md.columns:
                        if c.lower() in ('compound','pert_iname','perturbagen','sig_id'):
                            comp_col = c
                            break
                    val_col = None
                    for c in md.columns:
                        if c.lower() in ('value','zscore','score','logfc'):
                            val_col = c
                            break
                    if gene_col and comp_col and val_col:
                        wide = md.pivot_table(index=gene_col, columns=comp_col, values=val_col, aggfunc='mean')
                        wide.index = wide.index.astype(str).str.upper().str.strip()
                        lincs_df = wide
                        logging.info(f"Loaded manual_drug_scores and pivoted to wide (shape {wide.shape})")
                else:

                    try:
                        md2 = pd.read_csv(f, sep='\t', index_col=0)
                        md2.index = md2.index.astype(str).str.upper().str.strip()
                        lincs_df = md2
                        logging.info(f"Loaded manual_drug_scores as wide (shape {md2.shape})")
                    except Exception:
                        logging.error("Could not parse manual_drug_scores.tsv")
            except Exception as e:
                logging.error(f"Error reading manual_drug_scores.tsv: {e}")

    if lincs_df is None:
        logging.error("No LINCS data available. Provide --lincs PATH to a LINCS matrix.")
        return

    scores = compute_scores(sig, lincs_df)
    if scores is None:
        return

    scored = label_by_percentile(scores, percentile_cut=args.percentile)

    if args.smiles:
        sp = Path(args.smiles)
        if sp.exists():
            smdf = pd.read_csv(sp, sep='\t', index_col=0)

            joined = scored.join(smdf[['smiles']], how='left')
            scored = joined

    scored.to_csv(args.output, sep='\t')
    logging.info(f"Wrote labeled compounds to {args.output} (rows={len(scored)})")

if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--signature', required=True, default="C4A_complement_signature.tsv", help='Path to signature TSV (index=gene)')
    ap.add_argument('--lincs', default=r"C:\Users\pranj\Downloads\GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx\GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx", help='Path to LINCS L1000 wide matrix (genes x compounds) or long table')
    ap.add_argument('--output', default='labeled_compounds_C4A.tsv', help='Output labeled TSV')
    ap.add_argument('--percentile', type=float, default=5.0, help='Percentile cut for most-negative spearman to label reversers (default 5.0)')
    ap.add_argument('--topn', type=int, default=100, help='Number of genes to write to up_signature/down_signature')
    ap.add_argument('--smiles', help='Optional TSV with compound id as index and column smiles to attach to output')
    args = ap.parse_args()
    main(args)
