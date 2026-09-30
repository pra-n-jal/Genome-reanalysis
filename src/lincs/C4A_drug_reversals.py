from pathlib import Path
import argparse
import sys
import pandas as pd

def to_numeric_if_exists(df, col):
    if col in df.columns:
        try:
            df[col] = pd.to_numeric(df[col], errors='coerce')
        except Exception:
            pass

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--pert-scores', default='Reversal/connectivity_by_pert.tsv')
    parser.add_argument('--col-meta', default='lincs_col_meta_parsed.tsv')
    parser.add_argument('--pert-info', default='GSE92742_Broad_LINCS_pert_info.txt')
    parser.add_argument('--out', default='Reversal/C4A_drug_reversals.tsv')
    args = parser.parse_args()

    print('[DEBUG] Loading input files...')
    pert_scores = pd.read_csv(args.pert_scores, sep='\t', dtype=str, low_memory=False)
    col_meta = pd.read_csv(args.col_meta, sep='\t', dtype=str, low_memory=False)
    pert_info = pd.read_csv(args.pert_info, sep='\t', dtype=str, low_memory=False)
    to_numeric_if_exists(pert_scores, 'spearman')
    to_numeric_if_exists(pert_scores, 'reversal_score')
    print(f'[DEBUG] pert_scores columns: {list(pert_scores.columns)}')
    print(f'[DEBUG] col_meta columns:   {list(col_meta.columns)}')
    print(f'[DEBUG] pert_info columns:  {list(pert_info.columns)}')
    pert_scores.columns = [c.strip() for c in pert_scores.columns]
    col_meta.columns = [c.strip() for c in col_meta.columns]
    pert_info.columns = [c.strip() for c in pert_info.columns]
    if 'pert_id' not in pert_scores.columns:
        if 'sig_id' in pert_scores.columns:
            print("[DEBUG] 'pert_id' not found in pert_scores — attempting to map via col_meta using sig_id...")
            merged_tmp = pert_scores.merge(col_meta, left_on='sig_id', right_on='id', how='left')
            if 'full_pert_batch' in merged_tmp.columns and merged_tmp['full_pert_batch'].notna().any():
                merged_tmp['pert_id'] = merged_tmp['full_pert_batch']
                pert_scores = merged_tmp
                print('[DEBUG] Derived pert_id from full_pert_batch in col_meta.')
            elif 'pert_id' in merged_tmp.columns and merged_tmp['pert_id'].notna().any():
                pert_scores = merged_tmp
                print("[DEBUG] Derived pert_id from col_meta 'pert_id' column.")
            else:
                print('[FATAL] Could not derive pert_id from col_meta. Check your lincs_col_meta_parsed.tsv.')
                print('[DEBUG] Sample merged_tmp columns:', list(merged_tmp.columns))
                sys.exit(1)
        else:
            print("[FATAL] pert_scores lacks both 'pert_id' and 'sig_id'. Cannot proceed.")
            sys.exit(1)
    else:
        print("[DEBUG] pert_scores already contains 'pert_id' — will use it (with validation checks).")
    pert_scores['pert_id'] = pert_scores['pert_id'].astype(str).fillna('')
    short_ids = pert_scores['pert_id'].apply(lambda x: '-' not in x and x != '')
    if short_ids.any() and 'full_pert_batch' in col_meta.columns:
        if 'sig_id' in pert_scores.columns:
            print("[WARN] Detected short/ambiguous pert_id values (e.g., 'BRD'). Attempting to map to full_pert_batch via sig_id...")
            map_df = col_meta[['id', 'full_pert_batch']].dropna().drop_duplicates()
            map_dict = dict(zip(map_df['id'], map_df['full_pert_batch']))
            pert_scores['pert_id_full'] = pert_scores.apply(lambda r: map_dict.get(r.get('sig_id')) if pd.notna(r.get('sig_id')) else None, axis=1)
            replaced = pert_scores['pert_id_full'].notna().sum()
            if replaced > 0:
                pert_scores.loc[pert_scores['pert_id_full'].notna(), 'pert_id'] = pert_scores.loc[pert_scores['pert_id_full'].notna(), 'pert_id_full']
                print(f'[DEBUG] Replaced {replaced} short pert_id values with full_pert_batch.')
            else:
                print('[WARN] No sig_id -> full_pert_batch mappings found for short pert_id values. Proceeding, but merge may fail.')
            pert_scores.drop(columns=['pert_id_full'], inplace=True)
        else:
            print('[WARN] Short pert_id values detected but no sig_id column to map them. Proceeding (merge likely to fail).')
    has_brd = pert_info['pert_id'].astype(str).str.startswith('BRD').any()
    if not has_brd:
        print('[FATAL] Your pert_info appears not to contain any BRD- compound entries.')
        print('        It looks like a gene-only pert_info (trt_oe / trt_sh).')
        print('        Please supply the full GSE92742_Broad_LINCS_pert_info.txt that includes BRD- compound rows.')
        sys.exit(1)
    pert_info_cp = pert_info[pert_info['pert_type'] == 'trt_cp'].copy()
    print(f'[DEBUG] pert_info: total rows = {len(pert_info)}, compound (trt_cp) rows = {len(pert_info_cp)}')
    print('[DEBUG] Merging pert_scores with pert_info_cp on pert_id...')
    merged = pert_scores.merge(pert_info_cp, on='pert_id', how='left', indicator=True)
    print(f'[DEBUG] After merge: {merged.shape}')
    print(f"[DEBUG] Merge indicator counts:\n{merged['_merge'].value_counts().to_dict()}")
    left_only = merged[merged['_merge'] == 'left_only']
    if len(left_only) > 0:
        sample_left = left_only['pert_id'].dropna().unique()[:10].tolist()
        print(f'[WARN] {len(left_only)} rows did not match any compound in pert_info_cp. Sample unmatched pert_id(s): {sample_left}')
        if 'full_pert_batch' in col_meta.columns and 'sig_id' in merged.columns:
            print('[DEBUG] Attempting fallback mapping using col_meta.full_pert_batch where available...')
            map_df = col_meta[['id', 'full_pert_batch']].dropna().drop_duplicates()
            map_dict = dict(zip(map_df['id'], map_df['full_pert_batch']))
            fallback_mask = merged['_merge'] == 'left_only'
            count_fallback = 0
            for idx, row in merged[fallback_mask].iterrows():
                sid = row.get('sig_id')
                if pd.notna(sid) and sid in map_dict:
                    merged.at[idx, 'pert_id'] = map_dict[sid]
                    count_fallback += 1
            if count_fallback > 0:
                print(f'[DEBUG] Applied fallback mapping for {count_fallback} rows — re-merging...')
                print(f'[DEBUG] Applied fallback mapping for {count_fallback} rows — re-merging...')
                merged = merged.drop(columns=[c for c in merged.columns if c.endswith('_x') or c.endswith('_y') or c == '_merge'], errors='ignore')
                merged = merged.merge(pert_info_cp, on='pert_id', how='left', indicator=True)
                print(f'[DEBUG] After fallback re-merge: {merged.shape}')
                print(f"[DEBUG] Merge indicator counts after fallback: {merged['_merge'].value_counts().to_dict()}")
                print(f'[DEBUG] After fallback re-merge: {merged.shape}')
            else:
                print('[DEBUG] Fallback mapping found no replacements. If output lacks drugs, check your pert_info file.')
    if 'pert_iname' not in merged.columns or merged['pert_iname'].isna().all():
        print("[ERROR] No 'pert_iname' values found after merging with compound pert_info.")
        print('[DEBUG] Available columns after merge:', list(merged.columns))
        sys.exit(1)
    compound_data = merged[merged['pert_iname'].notna()].copy()
    if compound_data.empty:
        print('[ERROR] No compound perturbation rows with pert_iname present. Exiting.')
        sys.exit(1)
    print(f'[INFO] Retained {len(compound_data)} perturbation rows with compound names.')
    agg_cols = {}
    if 'spearman' in compound_data.columns:
        agg_cols['spearman'] = 'mean'
    if 'reversal_score' in compound_data.columns:
        agg_cols['reversal_score'] = 'mean'
    if 'single_gene_value' in compound_data.columns and compound_data['single_gene_value'].notna().any():
        agg_cols['single_gene_value'] = 'first'
    if not agg_cols:
        print('[WARN] No numeric columns (spearman/reversal_score/single_gene_value) found to aggregate — proceeding to produce pert_id/pert_iname counts only.')
    agg = compound_data.groupby(['pert_id', 'pert_iname']).agg(agg_cols if agg_cols else {'pert_id': 'size'}).reset_index()
    agg['n_signatures'] = compound_data.groupby(['pert_id', 'pert_iname']).size().values
    if 'reversal_score' in agg.columns:
        agg['reversal_score'] = pd.to_numeric(agg['reversal_score'], errors='coerce')
        agg = agg.sort_values('reversal_score', ascending=False)
    elif 'spearman' in agg.columns:
        agg['spearman'] = pd.to_numeric(agg['spearman'], errors='coerce')
        agg = agg.sort_values('spearman', ascending=False)
    outpath = Path(args.out)
    outpath.parent.mkdir(parents=True, exist_ok=True)
    agg.to_csv(outpath, sep='\t', index=False)
    print(f'[SUCCESS] Wrote pert-level (drug-level) ranked reversal table → {outpath}')
    print(f'[INFO] Total unique perturbagens ranked: {len(agg)}')
    print('\n=== Top 10 C4A Reversal Compounds (pert-level) ===')
    print(agg.head(10).to_string(index=False))


if __name__ == '__main__':
    main()
