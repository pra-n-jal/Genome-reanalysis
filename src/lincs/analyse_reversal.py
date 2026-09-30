import os
import sys
import argparse
import pandas as pd
import re
import math
from collections import Counter, defaultdict
try:
    from scipy.stats import hypergeom
    SCIPY_AVAILABLE = True
except Exception:
    SCIPY_AVAILABLE = False

def safe_read(path, **kwargs):
    try:
        return pd.read_csv(path, **kwargs)
    except Exception as e:
        print(f'[WARN] Could not read {path}: {e}')
        return None

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument('--project-root', '-p', default='.')
    p.add_argument('--reversal-tsv', default='Reversal/C4A_drug_reversals_with_pertid.tsv')
    p.add_argument('--pert-info', default='GSE92742_Broad_LINCS_pert_info.csv')
    p.add_argument('--expanded', default='Reversal/connectivity_by_pert_expanded.tsv')
    p.add_argument('--col-meta', default='lincs_col_meta_parsed.tsv')
    p.add_argument('--min-signatures', type=int, default=3)
    p.add_argument('--top-k', type=int, default=50)
    p.add_argument('--outdir', '-o', default='outputs')
    return p.parse_args()

def extract_pert_from_sig(sig_val):
    if pd.isna(sig_val):
        return None
    s = str(sig_val)
    tokens = s.split(':')
    for t in tokens:
        tt = t.strip()
        if re.search('\\b(BRD-|CMAP-|CMAP_|BRD_)\\b', tt, re.IGNORECASE) or tt.upper().startswith(('BRD', 'CMAP', 'CB', 'DRUG')):
            return tt
    m = re.search('(BRD-[A-Za-z0-9\\-]+)', s, re.IGNORECASE)
    if m:
        return m.group(1)
    m2 = re.search('(CMAP-[A-Za-z0-9\\-]+)', s, re.IGNORECASE)
    if m2:
        return m2.group(1)
    if len(tokens) >= 2:
        return tokens[1].strip()
    return s

def parse_dose_to_uM(dose_str):
    if pd.isna(dose_str):
        return None
    s = str(dose_str).replace(',', ' ').replace('μ', 'u').replace('µ', 'u').strip()
    m = re.search('([0-9.eE+-]+)\\s*([munp]?u?M|M|mM|nM|pM)?', s)
    if not m:
        try:
            return float(s)
        except:
            return None
    val = float(m.group(1))
    unit = (m.group(2) or '').lower()
    if unit in ('um', 'uM', 'µm', 'μm'):
        return float(val)
    if unit in ('nm', 'nm'):
        return float(val) / 1000.0
    if unit in ('mm', 'mm'):
        return float(val) * 1000.0
    if unit in ('m', 'm'):
        return float(val) * 1000000.0
    if unit in ('pm', 'pm'):
        return float(val) / 1000000.0
    return float(val)

def hypergeom_pval(k, M, n, N):
    if SCIPY_AVAILABLE:
        return float(hypergeom.sf(k - 1, M, n, N))
    denom = math.comb(M, N)
    top = 0.0
    for i in range(max(0, k), min(n, N) + 1):
        top += math.comb(n, i) * math.comb(M - n, N - i)
    return top / denom if denom > 0 else 1.0

def split_moa_field(moa_field):
    if pd.isna(moa_field):
        return []
    return [p.strip().lower() for p in re.split('\\||;|,', str(moa_field)) if p.strip()]

def main():
    args = parse_args()
    root = args.project_root
    os.makedirs(args.outdir, exist_ok=True)
    rev_path = os.path.join(root, args.reversal_tsv)
    print(f'[INFO] Loading reversal file: {rev_path}')
    rev = safe_read(rev_path, sep='\t', dtype=str)
    if rev is None:
        print('[ERROR] reversal TSV missing. Exiting.')
        sys.exit(1)
    sig_col = None
    for cand in ('sig_id', 'sig', 'signature_id', 'sigid'):
        if cand in rev.columns:
            sig_col = cand
            break
    if sig_col is None:
        for c in rev.columns:
            if 'sig' in c.lower() and 'id' in c.lower():
                sig_col = c
                break
    if sig_col is None:
        print(f'[WARN] No signature id column found. Columns: {list(rev.columns)}')
    else:
        print(f"[INFO] Using signature column '{sig_col}' to extract pert_id when needed.")
    if 'pert_id' not in rev.columns:
        rev['pert_id'] = None
    need_extract_mask = rev['pert_id'].isna() | (rev['pert_id'].astype(str).str.strip() == '')
    if sig_col is not None and need_extract_mask.any():
        print(f'[INFO] Extracting pert_id from {sig_col} for {need_extract_mask.sum()} rows')
        rev.loc[need_extract_mask, 'pert_id'] = rev.loc[need_extract_mask, sig_col].apply(extract_pert_from_sig)
    corrected_path = os.path.join(root, 'Reversal', 'C4A_drug_reversals_with_pertid.tsv')
    rev.to_csv(corrected_path, sep='\t', index=False)
    print(f'[INFO] Wrote corrected reversal TSV -> {corrected_path}')
    for col in ('reversal_score', 'single_gene_value', 'n_signatures'):
        if col in rev.columns:
            rev[col] = pd.to_numeric(rev[col], errors='coerce')
    pert_path = os.path.join(root, args.pert_info)
    pert = safe_read(pert_path, sep=None, engine='python', dtype=str)
    if pert is None:
        print(f'[WARN] Pert info {pert_path} not found or unreadable. Mapping will be limited.')
        pert = pd.DataFrame(columns=['pert_id', 'pert_iname', 'moa', 'pubchem_cid', 'drug_status', 'smiles'])
    pert_id_col = None
    for c in pert.columns:
        if c.lower() == 'pert_id' or 'pert_id' in c.lower():
            pert_id_col = c
            break
    if pert_id_col is None:
        for c in pert.columns:
            if c.lower() in ('pert_iname', 'pert_iname_label', 'perturbagen'):
                pert_id_col = c
                break
    if pert_id_col is None:
        print(f'[WARN] No pert id-like column found in pert info. Columns: {list(pert.columns)}')
    else:
        print(f"[INFO] Using pert-info column '{pert_id_col}' as join key.")
    if pert_id_col is not None:
        merged = rev.merge(pert, left_on='pert_id', right_on=pert_id_col, how='left', suffixes=('_rev', '_pert'))
    else:
        merged = rev.copy()
    merged['mapped'] = ~merged.filter(regex='moa|pubchem|drug_status|smiles', axis=1).isna().all(axis=1)
    if 'n_signatures' in merged.columns:
        merged['n_signatures'] = pd.to_numeric(merged['n_signatures'], errors='coerce').fillna(0).astype(int)
    else:
        merged['n_signatures'] = 0
    filtered = merged[merged['n_signatures'] >= args.min_signatures].copy()
    print(f'[INFO] Rows total: {len(merged)}; after n_signatures >= {args.min_signatures}: {len(filtered)}')
    expanded_path = os.path.join(root, args.expanded)
    expanded = safe_read(expanded_path, sep='\t', dtype=str)
    if expanded is not None:
        if 'pert_id' not in expanded.columns:
            for cand in ('sig_id', 'sig', 'signature_id'):
                if cand in expanded.columns:
                    expanded['pert_id'] = expanded[cand].apply(extract_pert_from_sig)
                    break
        if 'cell_id' in expanded.columns:
            cell_counts = expanded.groupby('pert_id')['cell_id'].nunique()
        else:
            cell_counts = pd.Series(dtype=int)
        sig_counts = expanded.groupby('pert_id').size()
        if 'dose' in expanded.columns:
            expanded['dose_uM'] = expanded['dose'].apply(parse_dose_to_uM)
            median_dose = expanded.groupby('pert_id')['dose_uM'].median()
        else:
            median_dose = pd.Series(dtype=float)
        filtered = filtered.set_index('pert_id')
        filtered['n_cell_lines'] = cell_counts.reindex(filtered.index).fillna(0).astype(int)
        filtered['n_expanded_signatures'] = sig_counts.reindex(filtered.index).fillna(0).astype(int)
        filtered['median_dose_uM'] = median_dose.reindex(filtered.index)
        filtered = filtered.reset_index()
    else:
        filtered['n_cell_lines'] = 0
        filtered['n_expanded_signatures'] = 0
        filtered['median_dose_uM'] = None
    filtered['reversal_score'] = pd.to_numeric(filtered.get('reversal_score', pd.Series([0] * len(filtered))), errors='coerce').fillna(0)
    filtered['reversal_rank'] = filtered['reversal_score'].rank(ascending=False, method='min')
    filtered['sig_rank'] = filtered['n_signatures'].rank(ascending=False, method='min')
    filtered['mapped_int'] = filtered['mapped'].astype(int)
    weights = {'rev': 0.6, 'sig': 0.3, 'mapped': 0.1}
    filtered['priority'] = filtered['reversal_rank'] * (1 - weights['rev']) + filtered['sig_rank'] * (1 - weights['sig']) - filtered['mapped_int'] * (weights['mapped'] * 10)
    filtered = filtered.sort_values('priority')
    topk = int(min(args.top_k, len(filtered)))
    top_df = filtered.head(topk).copy()
    all_df = filtered.copy()
    bg_moa_terms = Counter()
    top_moa_terms = Counter()
    for _, row in all_df.iterrows():
        for m in split_moa_field(row.get('moa', '')):
            bg_moa_terms[m] += 1
    for _, row in top_df.iterrows():
        for m in split_moa_field(row.get('moa', '')):
            top_moa_terms[m] += 1
    M = len(all_df)
    N = topk
    enrichment_records = []
    for moa, k in top_moa_terms.items():
        n = bg_moa_terms.get(moa, 0)
        pval = hypergeom_pval(int(k), M, int(n), N) if n > 0 else 1.0
        enrichment_records.append({'moa': moa, 'k_in_top': int(k), 'n_in_bg': int(n), 'M': int(M), 'N': int(N), 'pval': pval})
    moa_enrich_df = pd.DataFrame(enrichment_records).sort_values('pval')

    def is_high_dose(row):
        try:
            d = float(row.get('median_dose_uM'))
            return d > 10.0
        except:
            return False
    summary = filtered[['pert_id', 'pert_iname', 'reversal_score', 'single_gene_value', 'n_signatures', 'n_cell_lines', 'n_expanded_signatures', 'median_dose_uM', 'moa']].copy()
    summary['flag_high_dose'] = summary.apply(is_high_dose, axis=1)
    summary['flag_few_cell_lines'] = summary['n_cell_lines'].fillna(0).astype(int) <= 1
    if len(moa_enrich_df) > 0:
        top_moa = moa_enrich_df.iloc[0]['moa']
        summary['flag_moa_bias'] = summary['moa'].fillna('').str.lower().str.contains(str(top_moa))
    else:
        summary['flag_moa_bias'] = False
    out_mapped = os.path.join(args.outdir, 'mapped_prioritized.tsv')
    out_moa = os.path.join(args.outdir, 'moa_enrichment.tsv')
    out_summary = os.path.join(args.outdir, 'summary_flags.tsv')
    os.makedirs(args.outdir, exist_ok=True)
    filtered.to_csv(out_mapped, sep='\t', index=False)
    moa_enrich_df.to_csv(out_moa, sep='\t', index=False)
    summary.to_csv(out_summary, sep='\t', index=False)
    print(f'[DONE] wrote mapped/prioritized -> {out_mapped}')
    print(f'[DONE] wrote moa enrichment -> {out_moa}')
    print(f'[DONE] wrote summary flags -> {out_summary}')
    print('\nTop 20 prioritized:')
    show_cols = ['pert_id', 'pert_iname', 'reversal_score', 'n_signatures', 'n_cell_lines', 'median_dose_uM', 'moa']
    print(filtered[show_cols].head(20).to_string(index=False))
if __name__ == '__main__':
    main()
