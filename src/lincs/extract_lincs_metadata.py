import argparse
from pathlib import Path

import pandas as pd
from cmapPy.pandasGEXpress import parse_gctx


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--gctx', required=True, help='Path to the GCTX file')
    parser.add_argument('--out', default='lincs_col_meta.tsv')
    args = parser.parse_args()

    print('Reading GCTX column metadata...')
    gct = parse_gctx.parse(
        args.gctx,
        rid=None,
        cid=None,
        convert_neg_666=True,
        make_multiindex=False
    )

    if hasattr(gct, 'col_meta') and not gct.col_meta.empty:
        col_meta = gct.col_meta.copy()
    elif hasattr(gct, 'col_metadata_df') and not gct.col_metadata_df.empty:
        col_meta = gct.col_metadata_df.copy()
    else:
        raise ValueError('No column metadata found in GCTX. The file may be truncated or unreadable.')

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    col_meta.reset_index(drop=True).to_csv(out, sep='\t', index=False)

    print('Saved:', out)
    print('Columns:', list(col_meta.columns)[:12])
    print('Rows:', len(col_meta))


if __name__ == '__main__':
    main()
