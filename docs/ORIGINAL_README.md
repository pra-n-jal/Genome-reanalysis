# C4A / LINCS scripts

These are the cleaned versions of the scripts from the pipeline. The main calculations are left alone. The biggest changes are that the code is less cluttered, comments are kept short, and file locations are passed as arguments instead of being tied to one computer.

## Install

```bash
pip install pandas numpy scipy scikit-learn h5py requests tqdm
pip install mygene
```

For the LINCS GCTX scripts:

```bash
pip install cmapPy
```

For the chemistry scripts, install RDKit using the method that matches your Python setup. Conda is usually the least painful option:

```bash
conda install -c conda-forge rdkit
```

## Scripts

### 1. `gctx_tools.py`

Look at the GCTX structure:

```bash
python gctx_tools.py --mode explore --gctx /path/to/file.gctx
```

Extract column metadata:

```bash
python gctx_tools.py --mode extract \
    --gctx /path/to/file.gctx \
    --out lincs_col_meta.tsv
```

Parse the signature IDs from the extracted metadata:

```bash
python gctx_tools.py --mode parse \
    --input lincs_col_meta.tsv \
    --out lincs_col_meta.tsv
```

The parse step writes a file ending in `_parsed.tsv`.

### 2. `extract_lincs_metadata.py`

This is the shorter version that only extracts the GCTX column metadata.

```bash
python extract_lincs_metadata.py \
    --gctx /path/to/file.gctx \
    --out lincs_col_meta.tsv
```

### 3. `compute_connectivity_lincs.py`

Use a gene signature against LINCS. It can read a wide/long TSV or a GCTX file.

Using GCTX:

```bash
python compute_connectivity_lincs.py \
    --gctx /path/to/LINCS.gctx \
    --sig /path/to/signature.tsv \
    --pert_info /path/to/GSE92742_Broad_LINCS_pert_info.txt \
    --out Reversal/connectivity_by_pert.tsv
```

Using an already prepared LINCS table:

```bash
python compute_connectivity_lincs.py \
    --lincs /path/to/lincs.tsv \
    --sig /path/to/signature.tsv \
    --pert_info /path/to/GSE92742_Broad_LINCS_pert_info.txt \
    --out Reversal/connectivity_by_pert.tsv
```

The script also writes the per-signature table. Change `--out-sid` if you want it somewhere else.

### 4. `C4A_drug_reversals.py`

Join the connectivity results to LINCS perturbation information and produce the drug-level reversal table.

```bash
python C4A_drug_reversals.py \
    --pert-scores Reversal/connectivity_by_pert.tsv \
    --col-meta lincs_col_meta_parsed.tsv \
    --pert-info GSE92742_Broad_LINCS_pert_info.txt \
    --out Reversal/C4A_drug_reversals.tsv
```

### 5. `analyse_reversal.py`

This is the follow-up filtering, ranking and MOA analysis.

```bash
python analyse_reversal.py \
    --project-root /path/to/project \
    --reversal-tsv Reversal/C4A_drug_reversals_with_pertid.tsv \
    --pert-info GSE92742_Broad_LINCS_pert_info.csv \
    --expanded Reversal/connectivity_by_pert_expanded.tsv \
    --col-meta lincs_col_meta_parsed.tsv \
    --min-signatures 3 \
    --top-k 50 \
    --outdir outputs
```

### 6. `filter_and_cluster.py`

Calculate molecular descriptors, flag PAINS/nitro compounds and assign Butina clusters.

```bash
python filter_and_cluster.py \
    --input chemprop_C4A_predictions.csv \
    --out Reversal/predictions_filtered.tsv
```

The cluster cutoff can be changed with `--cluster-cutoff`.

### 7. `clean_predictions.py`

Apply the existing final filters:

- no PAINS flag
- no nitro group
- molecular weight below 600
- QED above 0.25

Run:

```bash
python clean_predictions.py \
    --input Reversal/predictions_filtered.tsv \
    --output Reversal/predictions_filtered_top.tsv
```

### 8. `map_to_lincs.py`

Map predicted molecules to the LINCS perturbation table by InChIKey.

```bash
python map_to_lincs.py \
    --pred chemprop_C4A_predictions.csv \
    --pert-info GSE92742_Broad_LINCS_pert_info.txt \
    --out Reversal/predictions_mapped_lincs.tsv
```

### 9. `prep_and_split_for_chemprop.py`

Prepare the training data and make the train/validation/test split.

```bash
python prep_and_split_for_chemprop.py \
    --input C4A_training_dataset_filled.tsv \
    --out-prefix chemprop \
    --seed 42
```

To upsample the positive class in the training set:

```bash
python prep_and_split_for_chemprop.py \
    --input C4A_training_dataset_filled.tsv \
    --out-prefix chemprop \
    --upsample \
    --seed 42
```

### 10. `reversal.py`

Combine Chemprop predictions with the C4A gene table and calculate the reversal scores.

```bash
python reversal.py \
    --pred chemprop_C4A_predictions.csv \
    --c4 Reversal/C4A_genes_combined.tsv \
    --out Reversal/reversal_ranked.tsv \
    --topk 200
```

### 11. `quick_sanity.py`

A quick check that the reversal output joins to the LINCS perturbation information.

```bash
python quick_sanity.py \
    --reversal Reversal/C4A_drug_reversals.tsv \
    --pert-info GSE92742_Broad_LINCS_pert_info.txt \
    --out Reversal/top20_annotated.tsv
```

## Typical order

A normal run is roughly:

```text
GCTX
  -> gctx_tools.py / extract_lincs_metadata.py
  -> compute_connectivity_lincs.py
  -> C4A_drug_reversals.py
  -> analyse_reversal.py
```

For the Chemprop side:

```text
C4A training data
  -> prep_and_split_for_chemprop.py
  -> Chemprop
  -> filter_and_cluster.py
  -> clean_predictions.py
  -> map_to_lincs.py
  -> reversal.py
```

The scripts are intentionally separate. That makes it easier to rerun one step without having to rerun the whole pipeline.

## Reproducibility

Use the same input files, Python environment and random seed when repeating the Chemprop data split. The scripts no longer assume a specific Windows username or Downloads folder. The paths in the examples are only examples; replace them with the actual files on your machine.
