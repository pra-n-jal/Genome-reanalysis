import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Descriptors, QED, rdMolDescriptors
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.ML.Cluster import Butina


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', default='chemprop_C4A_predictions.csv')
    parser.add_argument('--out', default='Reversal/predictions_filtered.tsv')
    parser.add_argument('--cluster-cutoff', type=float, default=0.2)
    args = parser.parse_args()

    pred = pd.read_csv(args.input).dropna(subset=['smiles'])

    params = FilterCatalogParams()
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
    catalog = FilterCatalog(params)
    nitro = Chem.MolFromSmarts('[N+](=O)[O-]')

    records = []
    fps = []

    for _, row in pred.iterrows():
        smiles = row['smiles']
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            continue

        scaffold = MurckoScaffold.GetScaffoldForMol(mol)
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, 2048)
        fps.append(fp)

        records.append({
            'smiles': smiles,
            'label': row['label'],
            'mw': Descriptors.MolWt(mol),
            'logp': Descriptors.MolLogP(mol),
            'hbd': rdMolDescriptors.CalcNumHBD(mol),
            'hba': rdMolDescriptors.CalcNumHBA(mol),
            'tpsa': rdMolDescriptors.CalcTPSA(mol),
            'qed': QED.qed(mol),
            'is_pains': len(list(catalog.GetMatches(mol))) > 0,
            'has_nitro': bool(mol.HasSubstructMatch(nitro)),
            'scaffold': Chem.MolToSmiles(scaffold) if scaffold else ''
        })

    df = pd.DataFrame(records)

    np_fps = []
    for fp in fps:
        arr = np.zeros((2048,), dtype=np.int32)
        from rdkit.DataStructs.cDataStructs import ConvertToNumpyArray
        ConvertToNumpyArray(fp, arr)
        np_fps.append(arr)

    rd_fps = [
        AllChem.GetMorganFingerprintAsBitVect(
            Chem.MolFromSmiles(smiles), 2, 2048
        )
        for smiles in df['smiles']
    ]

    dists = []
    for i in range(len(rd_fps)):
        for j in range(i + 1, len(rd_fps)):
            dists.append(1 - DataStructs.TanimotoSimilarity(rd_fps[i], rd_fps[j]))

    if len(rd_fps) > 0:
        clusters = Butina.ClusterData(
            dists,
            len(rd_fps),
            args.cluster_cutoff,
            isDistData=True
        )
    else:
        clusters = []

    cluster_ids = {}
    for cluster_id, cluster in enumerate(clusters):
        for index in cluster:
            cluster_ids[index] = cluster_id

    df['cluster_id'] = [cluster_ids.get(i, -1) for i in range(len(df))]

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, sep='\t', index=False)
    print('Wrote', out)


if __name__ == '__main__':
    main()
