# Antigravity UI / backend map

## Recommended screens

### 1. Dashboard
- Current run
- Last successful stage
- Input/output counts
- Error log

### 2. Data & Inputs
- GEO tables
- LINCS GCTX
- LINCS perturbation metadata
- Chemical metadata

### 3. Signature Builder
Backend:
- `src/signature/build_c_4_a_complement_signature.py`
- `src/signature/get_c4a_genes.py`
- `src/signature/convert_signature_to_entr.py`

Display:
- datasets used
- genes retained
- up/down signature
- missing mappings

### 4. LINCS Connectivity
Backend:
- `src/lincs/gctx_tools.py`
- `src/lincs/extract_lincs_metadata.py`
- `src/lincs/compute_connectivity_lincs.py`
- `src/lincs/C4A_drug_reversals.py`

Display:
- perturbation count
- connectivity distribution
- top/bottom compounds
- cell-line/dose metadata

### 5. Reversal Analysis
Backend:
- `src/lincs/analyse_reversal.py`
- `src/lincs/quick_sanity.py`
- `src/lincs/manual_drug_lookup.py`

Display:
- candidate table
- MOA analysis
- evidence/provenance

### 6. Chemprop
Backend:
- `src/chemprop/prepare_c4a_for_training.py`
- `src/chemprop/prep_and_split_for_chemprop.py`

Chemprop itself is an external training step.

Display:
- dataset size
- class balance
- split sizes
- model metrics
- prediction confidence

### 7. Chemical Filtering
Backend:
- `src/chemprop/filter_and_cluster.py`
- `src/chemprop/clean_predictions.py`
- `src/chemprop/fetch_smiles_pubchem.py`

Display:
- retained/rejected molecules
- rejection reason
- descriptors
- clusters

### 8. LINCS Mapping
Backend:
- `src/lincs/map_to_lincs.py`
- `src/lincs/map_lincs_to_smiles.py`
- `src/lincs/reversal.py`

Display:
- mapped compounds
- unmatched compounds
- reversal evidence

## UI principle

Keep the UI thin. Python scripts remain the computational backend. Every run should record parameters and generated files so that a result shown in the UI can be traced back to an input and command.
