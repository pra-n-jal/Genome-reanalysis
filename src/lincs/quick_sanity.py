import argparse
from pathlib import Path
import pandas as pd


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--reversal', default='Reversal/C4A_drug_reversals.tsv')
    parser.add_argument('--pert-info', default='GSE92742_Broad_LINCS_pert_info.txt')
    parser.add_argument('--out', default='Reversal/top20_annotated.tsv')
    args = parser.parse_args()

    rev = pd.read_csv(args.reversal, sep='\t')
    pert_info = pd.read_csv(args.pert_info, sep='\t', low_memory=False)

    print('Total rows in reversal:', len(rev))
    print('Unique pert_iname:', rev['pert_iname'].nunique())

    merged = rev.merge(
        pert_info[['pert_iname', 'pert_type', 'moa', 'pert_id']],
        on='pert_iname',
        how='left'
    )

    mapped = merged['pert_type'].notnull().sum()
    print(f'Mapped to pert_info: {mapped} / {len(merged)} ({mapped / len(merged):.1%})')
    print('pert_type counts:\n', merged['pert_type'].value_counts(dropna=False).head(20))

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    merged.sort_values('spearman').head(20).to_csv(out, sep='\t', index=False)
    print('Wrote', out)


if __name__ == '__main__':
    main()
