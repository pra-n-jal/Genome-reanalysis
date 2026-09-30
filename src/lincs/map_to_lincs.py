import argparse
from pathlib import Path

import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem


def canonical_inchikey(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None, None
        canonical = Chem.MolToSmiles(mol, isomericSmiles=True)
        inchikey = AllChem.MolToInchiKey(mol)
        return canonical, inchikey
    except Exception:
        return None, None


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--pred', default='chemprop_C4A_predictions.csv')
    parser.add_argument('--pert-info', default='GSE92742_Broad_LINCS_pert_info.txt')
    parser.add_argument('--out', default='Reversal/predictions_mapped_lincs.tsv')
    args = parser.parse_args()

    pred = pd.read_csv(args.pred)
    pred[['can_smiles', 'inchikey']] = pred['smiles'].apply(
        lambda s: pd.Series(canonical_inchikey(s))
    )
    pred = pred.dropna(subset=['inchikey'])

    pert = pd.read_csv(args.pert_info, sep='\t', low_memory=False)

    smiles_col = None
    for col in pert.columns:
        if 'smiles' in col.lower() or 'canonical' in col.lower():
            smiles_col = col
            break

    if smiles_col:
        pert['can_smiles'] = pert[smiles_col].astype(str)
        pert['inchikey'] = pert['can_smiles'].apply(
            lambda s: canonical_inchikey(s)[1] if pd.notna(s) and s else None
        )
    else:
        pert['inchikey'] = None

    merged = pred.merge(
        pert[['pert_id', 'pert_iname', 'inchikey']],
        on='inchikey',
        how='left'
    )

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(out, sep='\t', index=False)

    print('Wrote', out)
    print(merged[['smiles', 'label', 'pert_id', 'pert_iname', 'inchikey']].head(50))


if __name__ == '__main__':
    main()
