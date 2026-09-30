import argparse, os, re, difflib
from pathlib import Path
import pandas as pd, numpy as np
from scipy.stats import spearmanr
from cmapPy.pandasGEXpress import parse

def read_signature(sig_path):
    df = pd.read_csv(sig_path, sep=None, engine='python')
    cols = [c.strip().lower() for c in df.columns]
    gene_col = None
    val_col = None
    for c in df.columns:
        lc = c.strip().lower()
        if lc in ('gene', 'gene_symbol', 'symbol', 'gene symbol'):
            gene_col = c
        if lc in ('logfc', 'log_fc', 'combined_z', 'value', 'score', 'mean_logfc', 'mean_logFC'):
            val_col = c
    if gene_col is None:
        gene_col = df.columns[0]
    if val_col is None:
        if len(df.columns) >= 2:
            val_col = df.columns[1]
        else:
            raise ValueError('Signature must have two columns: gene_symbol and numeric value.')
    sig = df[[gene_col, val_col]].copy()
    sig.columns = ['gene_symbol', 'value']
    sig['gene_symbol'] = sig['gene_symbol'].astype(str).str.strip()
    sig = sig.dropna(subset=['gene_symbol'])
    sig['value'] = pd.to_numeric(sig['value'], errors='coerce').fillna(0.0)
    sig = sig.groupby('gene_symbol', as_index=False)['value'].median()
    return sig

def read_lincs_wide_or_long(path):
    p = Path(path)
    if not p.exists():
        return None
    try:
        df = pd.read_csv(p, sep='\t', index_col=0)
        if df.shape[0] >= 1 and df.index.astype(str).str.contains('\\w').any():
            df.index = df.index.astype(str).str.upper().str.strip()
            return df
    except Exception:
        pass
    try:
        df2 = pd.read_csv(p, sep='\t')
        cols = [c.lower() for c in df2.columns]
        gene_col = None
        comp_col = None
        val_col = None
        for c in df2.columns:
            cl = c.lower()
            if cl in ('gene', 'genes', 'gene_symbol', 'symbol'):
                gene_col = c
            if cl in ('compound', 'pert_iname', 'perturbagen', 'drug', 'sig_id'):
                if comp_col is None:
                    comp_col = c
            if cl in ('value', 'zscore', 'logfc', 'score', 'signature'):
                val_col = c
        if gene_col and comp_col and val_col:
            wide = df2.pivot_table(index=gene_col, columns=comp_col, values=val_col, aggfunc='median')
            wide.index = wide.index.astype(str).str.upper().str.strip()
            return wide
    except Exception:
        pass
    return None

def detect_gene_rows(sig_genes, row_meta, data_index):
    sig_up = [s.upper() for s in sig_genes]
    if row_meta is not None and row_meta.shape[1] > 0:
        gene_cols = [c for c in row_meta.columns if 'gene' in c.lower() or 'symbol' in c.lower() or 'pr_gene' in c.lower()]
        for c in gene_cols:
            rm = row_meta[c].astype(str).str.strip().str.upper()
            mapping = {}
            for i, v in enumerate(rm):
                mapping.setdefault(v, []).append(i)
            matched_positions = []
            matched_genes = []
            for g in sig_up:
                if g in mapping:
                    matched_positions.extend(mapping[g])
                    matched_genes.append(g)
            if matched_positions:
                return (c, sorted(set(matched_positions)), matched_genes)
    idx = [str(x).strip().upper() for x in data_index]
    mapping_idx = {v: i for i, v in enumerate(idx)}
    matched_positions = []
    matched_genes = []
    for g in sig_up:
        if g in mapping_idx:
            matched_positions.append(mapping_idx[g])
            matched_genes.append(g)
    if matched_positions:
        return ('index_exact', sorted(set(matched_positions)), matched_genes)
    for i, v in enumerate(idx):
        for g in sig_up:
            if g in v and i not in matched_positions:
                matched_positions.append(i)
                matched_genes.append(g)
    if matched_positions:
        return ('index_substring', sorted(set(matched_positions)), matched_genes)
    tokens = []
    for s in sig_genes:
        parts = re.split('[/\\\\\\|;:,]+|\\s+///\\s+|\\s+//\\s+|\\s+', str(s))
        tokens += [p.strip().upper() for p in parts if p.strip()]
    tokens = list(dict.fromkeys(tokens))
    for i, v in enumerate(idx):
        for t in tokens:
            if t in v and i not in matched_positions:
                matched_positions.append(i)
                matched_genes.append(t)
    if matched_positions:
        return ('token_substring', sorted(set(matched_positions)), matched_genes)
    norm_idx = [re.sub('[^A-Z0-9]', '', x) for x in idx]
    norm_sig = [re.sub('[^A-Z0-9]', '', s.upper()) for s in sig_genes]
    mapping_norm = {v: i for i, v in enumerate(norm_idx)}
    for sg, nsg in zip(sig_genes, norm_sig):
        if nsg in mapping_norm:
            matched_positions.append(mapping_norm[nsg])
            matched_genes.append(sg.upper())
    if matched_positions:
        return ('normalized_exact', sorted(set(matched_positions)), matched_genes)
    candidates = idx
    for g in sig_up:
        close = difflib.get_close_matches(g, candidates, n=5, cutoff=0.85)
        for c in close:
            i = candidates.index(c)
            if i not in matched_positions:
                matched_positions.append(i)
                matched_genes.append(g)
    if matched_positions:
        return ('fuzzy', sorted(set(matched_positions)), matched_genes)
    return (None, [], [])

def auto_map_columns_to_pertinfo(col_ids, pert_info_df):
    best_col = None
    best_count = 0
    col_set = set(map(str, col_ids))
    for c in pert_info_df.columns:
        overlap = col_set.intersection(set(pert_info_df[c].dropna().astype(str)))
        if len(overlap) > best_count:
            best_count = len(overlap)
            best_col = c
    return (best_col, best_count)

def compute_connectivity(sig_series, lincs_wide, row_positions):
    genes = [lincs_wide.index[pos] for pos in row_positions]
    sig_map = {g.upper(): v for g, v in sig_series.items()}
    disease_values = []
    for g in genes:
        gu = str(g).upper()
        key = None
        for sg in sig_map:
            if sg == gu or sg in gu or gu in sg:
                key = sg
                break
        if key is None:
            key = list(sig_map.keys())[0]
        disease_values.append(sig_map[key])
    disease_values = np.array(disease_values, dtype=float)
    use_single = len(disease_values) < 3
    results = []
    for cid in lincs_wide.columns:
        try:
            if use_single:
                val = lincs_wide.iloc[row_positions[0], lincs_wide.columns.get_loc(cid)]
                rev = -float(val) if pd.notna(val) else np.nan
                results.append((cid, np.nan, float(val) if pd.notna(val) else np.nan, rev))
            else:
                expr = lincs_wide.iloc[row_positions, lincs_wide.columns.get_loc(cid)].values.astype(float)
                mask = ~np.isnan(expr) & ~np.isnan(disease_values)
                if mask.sum() < max(3, int(len(disease_values) * 0.1)):
                    results.append((cid, np.nan, np.nan, np.nan))
                    continue
                rho, _ = spearmanr(disease_values[mask], expr[mask])
                results.append((cid, float(rho) if not np.isnan(rho) else np.nan, np.nan, -float(rho) if not np.isnan(rho) else np.nan))
        except Exception:
            results.append((cid, np.nan, np.nan, np.nan))
    resdf = pd.DataFrame(results, columns=['sig_id', 'spearman', 'single_gene_value', 'reversal_score']).set_index('sig_id')
    return resdf

def main(args):
    print('[INFO] Loading signature:', args.sig)
    sig = read_signature(args.sig)
    print('[INFO] Signature genes:', len(sig))
    lincs_wide = None
    if args.lincs:
        lincs_wide = read_lincs_wide_or_long(args.lincs)
        if lincs_wide is not None:
            print('[INFO] Read LINCS via read_lincs_wide_or_long with shape', lincs_wide.shape)
    if lincs_wide is None and args.gctx:
        print('[INFO] Parsing GCTX via cmapPy:', args.gctx)
        gct = parse.parse(args.gctx)
        data_df = gct.data_df
        row_meta = gct.row_metadata_df
        col_meta = gct.col_metadata_df
        print('[INFO] GCTX loaded: data shape', data_df.shape, '; row_meta shape', getattr(row_meta, 'shape', None), '; col_meta shape', getattr(col_meta, 'shape', None))
        gene_col = None
        if row_meta is not None and row_meta.shape[1] > 0:
            for c in row_meta.columns:
                if 'gene' in c.lower() or 'symbol' in c.lower() or 'pr_gene' in c.lower():
                    gene_col = c
                    break
        if gene_col:
            new_index = row_meta[gene_col].astype(str).str.strip().str.upper().tolist()
            data_df.index = new_index
        else:
            data_df.index = data_df.index.astype(str)
        lincs_wide = data_df.copy()
        try:
            col_meta_df = col_meta
        except Exception:
            col_meta_df = pd.DataFrame(index=lincs_wide.columns)
    if lincs_wide is None:
        raise ValueError("Couldn't read LINCS data. Provide --lincs (wide or long) or --gctx.")
    lincs_wide.index = lincs_wide.index.astype(str).str.upper().str.strip()
    method, matched_positions, matched_genes = detect_gene_rows(sig['gene_symbol'].tolist(), gct.row_metadata_df if 'gct' in locals() else None, lincs_wide.index)
    if method is None or len(matched_positions) == 0:
        preview_n = min(500, lincs_wide.shape[0])
        preview_df = pd.DataFrame({'row_index': list(lincs_wide.index[:preview_n])})
        preview_out = Path(args.preview_out)
        preview_df.to_csv(preview_out, sep='\t', index=False)
        print('[ERROR] Could not match any signature genes to LINCS rows. Wrote preview:', preview_out)
        return
    print(f"[INFO] matched method '{method}' matched genes: {matched_genes} (n={len(matched_positions)})")
    resdf = compute_connectivity(sig.set_index('gene_symbol')['value'], lincs_wide, matched_positions)
    pert_map = None
    if args.pert_info and Path(args.pert_info).exists():
        pinfo = pd.read_csv(args.pert_info, sep='\t', dtype=str, low_memory=False)
        best_col, best_count = auto_map_columns_to_pertinfo(resdf.index.tolist(), pinfo)
        if best_col and best_count > 0:
            print(f"[INFO] Using pert_info column '{best_col}' (overlap={best_count}) to map signatures -> pert_info")
            pinfo = pinfo.set_index(best_col)
            possible_pid = [c for c in pinfo.columns if 'pert_id' in c.lower() or c.lower() == 'pertid']
            possible_iname = [c for c in pinfo.columns if 'pert_iname' in c.lower() or 'iname' in c.lower() or 'pert_iname' == c.lower()]
            pid_col = possible_pid[0] if possible_pid else None
            piname_col = possible_iname[0] if possible_iname else None
            pert_ids = []
            pert_inames = []
            for sid in resdf.index.astype(str):
                if sid in pinfo.index:
                    row = pinfo.loc[sid]
                    if isinstance(row, pd.DataFrame):
                        row = row.iloc[0]
                    pert_ids.append(row[pid_col] if pid_col and pid_col in row.index else '')
                    pert_inames.append(row[piname_col] if piname_col and piname_col in row.index else '')
                else:
                    pert_ids.append('')
                    pert_inames.append('')
            resdf['pert_id'] = pert_ids
            resdf['pert_iname'] = pert_inames
        else:
            print('[WARNING] Could not find an overlapping signature id column in pert_info; per-sid results will be written but cannot be aggregated by pert_id.')
            resdf['pert_id'] = ''
            resdf['pert_iname'] = ''
    else:
        resdf['pert_id'] = ''
        resdf['pert_iname'] = ''
    out_sid = Path(args.out_sid)
    resdf.to_csv(out_sid, sep='\t', index=True)
    out_sid.parent.mkdir(parents=True, exist_ok=True)
    print('[INFO] Wrote per-signature results to', out_sid)
    if resdf['pert_id'].astype(bool).any():
        agg = resdf.groupby('pert_id').agg(n_signatures=('reversal_score', 'count'), reversal_median=('reversal_score', 'median'), reversal_mean=('reversal_score', 'mean'), spearman_median=('spearman', 'median'), spearman_mean=('spearman', 'mean')).reset_index()
        agg.to_csv(args.out, sep='\t', index=False)
        print('[INFO] Wrote aggregated per-pert results to', args.out)
    elif resdf['pert_iname'].astype(bool).any():
        agg = resdf.groupby('pert_iname').agg(n_signatures=('reversal_score', 'count'), reversal_median=('reversal_score', 'median'), reversal_mean=('reversal_score', 'mean'), spearman_median=('spearman', 'median'), spearman_mean=('spearman', 'mean')).reset_index().rename(columns={'pert_iname': 'pert_iname'})
        agg.to_csv(args.out, sep='\t', index=False)
        print('[INFO] Wrote aggregated per-pert_iname results to', args.out)
    else:
        resdf.reset_index().to_csv(args.out, sep='\t', index=False)
        print('[INFO] No pert mapping available; wrote per-signature table to', args.out)
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--gctx', default=None, help='Path to .gctx (level5) file')
    parser.add_argument('--lincs', default=None, help='Optional wide or long LINCS file (wide gene x sig or long table)')
    parser.add_argument('--sig', required=True, help='Signature TSV (gene,value)')
    parser.add_argument('--shortlist', default=None, help='Optional shortlist (not required if pert_info mapping works)')
    parser.add_argument('--pert_info', default=None, help='GSE92742_Broad_LINCS_pert_info.txt (tab)')
    parser.add_argument('--out', default='Reversal/connectivity_by_pert.tsv', help='Output aggregated TSV')
    parser.add_argument('--out-sid', default='Reversal/connectivity_by_sid.tsv')
    parser.add_argument('--preview-out', default='Reversal/gctx_index_preview.tsv')
    parser.add_argument('--min_genes', type=int, default=3)
    args = parser.parse_args()
    main(args)
