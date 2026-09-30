import argparse

import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.utils import resample


def load_and_clean(path):
    df = pd.read_csv(path, sep='\t', low_memory=False)

    if 'smiles' not in df.columns:
        raise ValueError("Input missing 'smiles' column")
    if 'label_binary' not in df.columns:
        raise ValueError("Input missing 'label_binary' column. Run previous mapping step.")

    extra = [
        col for col in df.columns
        if col.startswith('spearman') or col.startswith('cosine')
    ]
    df = df[['smiles', 'score', 'label_binary'] + extra]

    df['smiles'] = df['smiles'].astype(str).str.strip()
    df = df[df['smiles'].notna() & (df['smiles'] != '')].copy()

    df = (
        df.groupby('smiles')
        .agg({'score': 'mean', 'label_binary': 'max'})
        .reset_index()
    )
    df['label_binary'] = df['label_binary'].astype(int)
    return df


def split_save(df, out_prefix, upsample=False, random_state=42):
    X = df['smiles']
    y = df['label_binary']

    X_train, X_temp, y_train, y_temp = train_test_split(
        X, y, test_size=0.2, stratify=y, random_state=random_state
    )
    X_val, X_test, y_val, y_test = train_test_split(
        X_temp, y_temp, test_size=0.5, stratify=y_temp, random_state=random_state
    )

    train_df = pd.DataFrame({'smiles': X_train, 'label': y_train}).reset_index(drop=True)
    val_df = pd.DataFrame({'smiles': X_val, 'label': y_val}).reset_index(drop=True)
    test_df = pd.DataFrame({'smiles': X_test, 'label': y_test}).reset_index(drop=True)

    if upsample:
        pos = train_df[train_df['label'] == 1]
        neg = train_df[train_df['label'] == 0]

        if len(pos) == 0:
            print('No positives in train to upsample.')
        elif len(pos) < len(neg):
            pos_up = resample(
                pos,
                replace=True,
                n_samples=len(neg),
                random_state=random_state
            )
            train_df = (
                pd.concat([neg, pos_up])
                .sample(frac=1, random_state=random_state)
                .reset_index(drop=True)
            )
            print(f'Upsampled positives: {len(pos)} -> {len(pos_up)}')

    train_df.to_csv(out_prefix + '_train.csv', index=False)
    val_df.to_csv(out_prefix + '_val.csv', index=False)
    test_df.to_csv(out_prefix + '_test.csv', index=False)

    full = df.rename(columns={'label_binary': 'label'}).copy()
    full.to_csv(out_prefix + '_full.csv', index=False)
    print('Saved train/val/test and full files.')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', default='C4A_training_dataset_filled.tsv')
    parser.add_argument('--out-prefix', default='chemprop')
    parser.add_argument('--upsample', action='store_true')
    parser.add_argument('--seed', type=int, default=42)
    args = parser.parse_args()

    df = load_and_clean(args.input)
    split_save(
        df,
        args.out_prefix,
        upsample=args.upsample,
        random_state=args.seed
    )


if __name__ == '__main__':
    main()
