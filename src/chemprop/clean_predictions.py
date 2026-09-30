import argparse
from pathlib import Path
import pandas as pd


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', default='Reversal/predictions_filtered.tsv')
    parser.add_argument('--output', default='Reversal/predictions_filtered_top.tsv')
    args = parser.parse_args()

    df = pd.read_csv(args.input, sep='\t')
    keep = df[
        (df['is_pains'] == False)
        & (df['has_nitro'] == False)
        & (df['mw'] < 600)
        & (df['qed'] > 0.25)
    ]

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    keep.sort_values('label', ascending=False).to_csv(out, sep='\t', index=False)
    print('Kept:', len(keep))
    print('Wrote', out)


if __name__ == '__main__':
    main()
