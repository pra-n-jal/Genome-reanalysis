
import argparse
import logging
import sys
from pathlib import Path
from collections import Counter, defaultdict

import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from tqdm import tqdm

logging.basicConfig(level=logging.INFO, format='[%(levelname)s] %(message)s')

def canonicalize_smiles(smi):

    if smi is None:
        return None, None
    s = str(smi).strip()
    if s == '' or s.lower() in ('nan', 'none', '-666'):
        return None, None
    try:
        m = Chem.MolFromSmiles(s, sanitize=True)
        if m is None:
            return None, None

        frags = Chem.GetMolFrags(m, asMols=True, sanitizeFrags=False)
        if len(frags) > 1:
            frags = sorted(frags, key=lambda x: x.GetNumHeavyAtoms(), reverse=True)
            m = frags[0]

        try:
            Chem.SanitizeMol(m)
        except Exception:
            pass
        can = Chem.MolToSmiles(m, isomericSmiles=True)
        try:
            ik = AllChem.MolToInchiKey(m)
        except Exception:
            ik = None
        return can, ik
    except Exception:
        return None, None

def build_report_and_clean(df, smiles_col, target_col, conflict_policy='majority', majority_threshold=0.7):

    if smiles_col not in df.columns:
        raise ValueError(f"SMILES column '{smiles_col}' not found in input.")
    if target_col not in df.columns:
        raise ValueError(f"Target column '{target_col}' not found in input.")

    tqdm.pandas(desc='Canonicalizing SMILES')
    results = df[smiles_col].progress_apply(canonicalize_smiles)
    df[['canonical_smiles', 'inchikey']] = pd.DataFrame(results.tolist(), index=df.index)

    mapped_cols = []
    for c in ['fetched_smiles_offline', 'fetched_smiles_pubchem', 'smiles', 'SMILES']:
        if c in df.columns:
            mapped_cols.append(c)

    def provenance_row(row):
        sources = set()
        if 'fetched_smiles_offline' in row.index and pd.notna(row.get('fetched_smiles_offline')) and str(row.get('fetched_smiles_offline')).strip():
            sources.add('repurposing')
        if 'fetched_smiles_pubchem' in row.index and pd.notna(row.get('fetched_smiles_pubchem')) and str(row.get('fetched_smiles_pubchem')).strip():
            sources.add('pubchem')

        if pd.notna(row.get(smiles_col)) and str(row.get(smiles_col)).strip():
            sources.add('original')
        if not sources:
            sources.add('unknown')
        return ';'.join(sorted(sources))

    df['mapped_source'] = df.apply(provenance_row, axis=1)

    before = len(df)
    df_valid = df[df['canonical_smiles'].notna()].copy()
    after = len(df_valid)
    logging.info(f"Dropped {before-after} rows with invalid SMILES (kept {after} rows)")

    key_col = 'inchikey'
    if df_valid['inchikey'].isna().all():
        key_col = 'canonical_smiles'
        logging.info('No InChIKeys produced; grouping by canonical_smiles instead of inchikey')
    else:
        logging.info('Grouping by InChIKey (preferred)')

    grouped = df_valid.groupby(key_col, dropna=False)

    clean_rows = []
    report_rows = []

    for key, group in tqdm(grouped, desc='Deduplicating groups'):
        canonical = group['canonical_smiles'].iloc[0]
        inchikey = group['inchikey'].iloc[0] if 'inchikey' in group.columns else None
        n = len(group)

        labels = group[target_col].dropna().astype(str).tolist()
        label_counts = Counter(labels)
        n_unique_labels = len(label_counts)
        sources = ';'.join(sorted(group['mapped_source'].dropna().unique().tolist()))

        decided_label = None
        kept = False
        if n_unique_labels == 0:

            report_rows.append((canonical, inchikey, n, n_unique_labels, '', sources, 'no_label', 'drop'))
            continue
        if n_unique_labels == 1:
            decided_label = labels[0]
            kept = True
            reason = 'consistent'
        else:
            if conflict_policy == 'majority':
                most_common, count = label_counts.most_common(1)[0]
                prop = count / sum(label_counts.values())
                if prop >= majority_threshold:
                    decided_label = most_common
                    kept = True
                    reason = f'majority({prop:.2f})'
                else:
                    decided_label = None
                    kept = False
                    reason = f'no_majority({prop:.2f})'
            elif conflict_policy == 'drop':
                decided_label = None
                kept = False
                reason = 'drop_conflict'
            elif conflict_policy == 'keep_first':
                decided_label = labels[0]
                kept = True
                reason = 'keep_first'
            else:
                decided_label = None
                kept = False
                reason = 'unknown_policy'

        report_rows.append((canonical, inchikey, n, n_unique_labels, dict(label_counts), sources, 'kept' if kept else 'dropped', reason))

        if kept:
            clean_rows.append({'canonical_smiles': canonical,
                               'smiles': canonical,
                               'inchikey': inchikey,
                               'target': decided_label,
                               'n_merged': n,
                               'merged_sources': sources})

    report_df = pd.DataFrame(report_rows, columns=['canonical_smiles', 'inchikey', 'n_rows', 'n_unique_labels', 'label_counts', 'merged_sources', 'status', 'reason'])
    clean_df = pd.DataFrame(clean_rows)

    return clean_df, report_df

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--input', '-i', default='C4A_training_dataset_filled.tsv', help='Input filled dataset (tab-separated)')
    ap.add_argument('--output', '-o', default='C4A_training_dataset_ready.tsv', help='Output cleaned dataset (tab-separated)')
    ap.add_argument('--report', default='dedup_report.csv', help='Dedup/report CSV path')
    ap.add_argument('--smiles-col', default='smiles', help='Name of the SMILES column in input')
    ap.add_argument('--target-col', required=True, help='Name of the target/label column in input')
    ap.add_argument('--conflict-policy', choices=['majority', 'drop', 'keep_first'], default='majority', help='Policy to resolve conflicting labels for same compound')
    ap.add_argument('--majority-threshold', type=float, default=0.7, help='Threshold for majority policy (0-1)')
    ap.add_argument('--min-rows', type=int, default=1, help='Minimum rows merged to keep group (for filtering tiny groups)')
    args = ap.parse_args()

    inp = Path(args.input)
    if not inp.exists():
        logging.error('Input file not found: %s', inp)
        sys.exit(1)

    logging.info('Loading input: %s', inp)
    df = pd.read_csv(inp, sep='\t', low_memory=False)
    logging.info('Loaded %d rows, columns: %s', len(df), ', '.join(df.columns.tolist()))

    clean_df, report_df = build_report_and_clean(df, args.smiles_col, args.target_col, conflict_policy=args.conflict_policy, majority_threshold=args.majority_threshold)

    out = Path(args.output)
    rep = Path(args.report)
    logging.info('Saving cleaned dataset to %s (rows=%d)', out, len(clean_df))

    cols = ['smiles', 'target', 'inchikey', 'n_merged', 'merged_sources']
    for c in cols:
        if c not in clean_df.columns:
            clean_df[c] = None
    clean_df[cols].to_csv(out, sep='\t', index=False)

    logging.info('Saving dedup report to %s (groups=%d)', rep, len(report_df))
    report_df.to_csv(rep, index=False)

    kept = len(clean_df)
    total_groups = len(report_df)
    dropped = report_df[report_df['status'] == 'dropped'].shape[0]
    logging.info('Summary: groups=%d, kept=%d, dropped=%d', total_groups, kept, dropped)
    logging.info('Done.')

if __name__ == '__main__':
    main()
