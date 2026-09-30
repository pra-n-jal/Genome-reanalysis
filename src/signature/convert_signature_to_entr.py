
import pandas as pd
from mygene import MyGeneInfo

sig = pd.read_csv("C4A_complement_signature.tsv", sep="\t")

mg = MyGeneInfo()
out = mg.querymany(sig['gene'].tolist(), scopes='symbol', fields='entrezgene', species='human')

mapping = {d['query']: d.get('entrezgene', None) for d in out}

sig['entrez_id'] = sig['gene'].map(mapping)
sig = sig.dropna(subset=['entrez_id'])

sig[['entrez_id', 'combined_z']].rename(columns={'entrez_id': 'gene'}).to_csv("C4A_signature_entrez.tsv", sep="\t", index=False)

print(f"[INFO] Saved mapped signature with {len(sig)} genes.")
