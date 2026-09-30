import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def load_predictions(path):
    df = pd.read_csv(path, sep=None, engine='python')
    df.columns = [c.strip() for c in df.columns]

    if 'label' not in df.columns:
        if df.shape[1] >= 2:
            df.columns = ['smiles', 'label'] + list(df.columns[2:])
        else:
            raise ValueError('Predictions file must contain SMILES and a label column.')

    df = df[['smiles', 'label']].copy()
    df['label'] = pd.to_numeric(df['label'], errors='coerce')
    return df.dropna(subset=['smiles'])


def load_c4_rows(path):
    df = pd.read_csv(path, sep=None, engine='python')
    columns = [c.lower() for c in df.columns]

    if 'dataset' not in columns:
        if 'gse' in df.columns[0].lower():
            df = df.rename(columns={df.columns[0]: 'dataset'})

    logfc_col = None
    for col in df.columns:
        name = col.lower()
        if name in ('logfc', 'log_fc', 'log2foldchange', 'log2fold_change'):
            logfc_col = col
            break

    if logfc_col is None:
        for col in df.columns:
            name = col.lower()
            if 'log' in name and 'fc' in name:
                logfc_col = col
                break

    if logfc_col is None:
        raise ValueError('Could not find a logFC column in the C4 file. Rename your logFC column to logFC.')

    if 'dataset' not in [c.lower() for c in df.columns]:
        raise ValueError('C4 file must contain a dataset column with GSE identifiers.')

    df = df.rename(columns={logfc_col: 'logFC'})
    df['logFC'] = pd.to_numeric(df['logFC'], errors='coerce').fillna(0.0)
    df['abs_logFC'] = df['logFC'].abs()
    return df[['dataset', 'gene_symbol', 'logFC', 'abs_logFC']]


def compute_scores(pred_df, c4_df):
    results = pred_df.groupby('smiles', as_index=False)['label'].mean()

    for _, row in c4_df.iterrows():
        score_col = f"score_{row['dataset']}"
        results[score_col] = results['label'] * row['abs_logFC']

    score_cols = [c for c in results.columns if c.startswith('score_')]
    score_values = results[score_cols].replace([np.inf, -np.inf], np.nan)

    results['score_mean'] = score_values.mean(axis=1)
    results['score_median'] = score_values.median(axis=1)
    results['score_sum'] = score_values.sum(axis=1)
    results['datasets_count'] = score_values.notna().sum(axis=1)

    results = results.sort_values(
        ['score_mean', 'score_sum'],
        ascending=False
    ).reset_index(drop=True)
    results['rank'] = results['score_mean'].rank(
        method='min', ascending=False
    ).astype(int)
    return results


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--pred', required=True)
    parser.add_argument('--c4', required=True)
    parser.add_argument('--out', default='reversal_ranked.tsv')
    parser.add_argument('--topk', type=int, default=200)
    args = parser.parse_args()

    pred = load_predictions(args.pred)
    c4 = load_c4_rows(args.c4)
    scored = compute_scores(pred, c4)

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    scored.to_csv(out, sep='\t', index=False)
    print(f'Wrote full ranked file: {out} (n={len(scored)})')

    if args.topk > 0:
        top_file = out.with_name(out.stem + f'_top{args.topk}.tsv')
        scored.head(args.topk).to_csv(top_file, sep='\t', index=False)
        print(f'Wrote top-{args.topk} file: {top_file}')


if __name__ == '__main__':
    main()
