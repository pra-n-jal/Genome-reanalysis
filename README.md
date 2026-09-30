# C4A / LINCS Research Pipeline

Welcome to the C4A / LINCS Research Project! This repository contains a complete bioinformatics and machine learning pipeline (using Chemprop) for analyzing C4A complement signatures and mapping them to drug reversals using the LINCS dataset.

## Getting Started

The easiest way to interact with the pipeline is via the included Desktop GUI.

### Prerequisites

Ensure you have Python 3.8+ installed, then install the required packages:

```bash
pip install -r requirements.txt
```

### Running the GUI

We've built a user-friendly, native desktop application to manage the pipeline scripts without needing to memorize command-line flags. 

```bash
python gui.py
```

* **Browse** for your input datasets (which should be placed in `data/raw/` or `data/intermediate/`).
* **Select** the script you want to run from the dropdown menu.
* Add any extra arguments (like `--out outputs/my_result.tsv`).
* **Run** the script and watch the logs in the integrated console!
* *New to the project?* Click **LOAD TEST EXAMPLE** in the GUI to see exactly how a run should be configured.

## Project Structure

* `gui.py` - The main desktop interface.
* `src/` - The core backend Python scripts (Signature Builder, LINCS Connectivity, Chemprop utilities).
* `data/` - Contains `raw/`, `intermediate/`, and `outputs/` subdirectories. *(Note: Large data files are ignored by git to keep the repository lightweight. You must provide your own GEO/LINCS datasets).*
* `docs/` - Contains the original README and in-depth scientific analysis.

For a detailed breakdown of the internal script workflow and architecture, see [PROJECT_STRUCTURE.md](PROJECT_STRUCTURE.md).
