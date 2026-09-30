# C4A / LINCS Research Project

This directory is an organized copy of the supplied research scripts and documentation.
The Python source files were copied without changing their contents.

## Directory layout

```text
C4A_LINCS_Project/
├── config/                  # Put run configuration here
├── data/
│   ├── raw/                 # GEO/LINCS/raw external inputs
│   ├── intermediate/        # Generated tables used between pipeline stages
│   └── outputs/             # Final tables, rankings and reports
├── docs/
│   ├── Research_In_Depth_Analysis.pdf
│   └── ORIGINAL_README.md
├── src/
│   ├── signature/           # C4A/complement signature construction + ID conversion
│   ├── lincs/               # LINCS connectivity, perturbation mapping and reversal analysis
│   ├── chemprop/            # Chemical training data, Chemprop prep and chemical filtering
│   └── utilities/           # Gene mapping and sanity checks
├── requirements.txt
└── PROJECT_STRUCTURE.md
```

## Intended workflow

### Signature branch

```text
GEO-derived gene tables
    ↓
src/signature/build_c_4_a_complement_signature.py
    ↓
C4A/complement signature
    ↓
src/signature/convert_signature_to_entr.py   (if Entrez IDs are required)
```

### Direct LINCS branch

```text
LINCS GCTX
    ↓
src/lincs/gctx_tools.py
    or
src/lincs/extract_lincs_metadata.py
    ↓
src/lincs/compute_connectivity_lincs.py
    ↓
src/lincs/C4A_drug_reversals.py
    ↓
src/lincs/analyse_reversal.py
```

### Chemprop branch

```text
C4A/LINCS-labelled training data
    ↓
src/chemprop/prepare_c4a_for_training.py
    ↓
src/chemprop/prep_and_split_for_chemprop.py
    ↓
Chemprop training (external)
    ↓
Chemprop predictions
    ↓
src/chemprop/filter_and_cluster.py
    ↓
src/chemprop/clean_predictions.py
    ↓
src/lincs/map_to_lincs.py
    ↓
src/lincs/reversal.py
```

### Supporting utilities

```text
src/utilities/map_genes.py
src/utilities/check_genes.py
src/chemprop/fetch_smiles_pubchem.py
src/lincs/map_lincs_to_smiles.py
src/lincs/quick_sanity.py
src/lincs/manual_drug_lookup.py
```

## UI integration notes

For a future Antigravity UI, treat the Python files as backend pipeline stages rather than as pages/components.
A useful UI model is:

1. **Project / Inputs** — select raw GEO/LINCS files.
2. **Signature** — build and inspect the C4A-associated signature.
3. **LINCS Connectivity** — run/query connectivity and inspect perturbations.
4. **Drug Reversal Analysis** — aggregate, filter and rank perturbagens.
5. **Chemprop** — prepare data, launch an external Chemprop run, inspect predictions.
6. **Chemical Filtering** — PAINS/QED/MW/cluster filtering.
7. **Candidate Mapping** — map predicted molecules back to LINCS.
8. **Results** — tables, plots, provenance and exported files.
9. **Run Log** — command, parameters, seed, input hashes and output locations.

The UI should invoke scripts through their CLI interfaces and display their generated files; it should not duplicate their scientific logic.

## Important distinction

The supplied research PDF describes a C4A-expression/correlation → LINCS connectivity workflow, while the supplied script assortment contains a later/alternate complement-signature + LINCS + Chemprop workflow. Keep these as separate documented workflow versions until the scientific methodology is reconciled.
