# BactoCat 

A computational framework to bridge the gap between *in vitro* and *in vivo* enzyme kinetics. **BactoCat** is an open-source, automated pipeline for computing apparent in vivo turnover numbers ($k_\text{app}$).

Given a genome-scale metabolic model (GEM), experimental conditions (uptake fluxes), and/or protein abundance data, it:
- **Parses GPR rules** — classifies every enzyme in the GEM as homomeric, 
  isoenzyme, or complex, and retains only those tractable for $k_\text{app}$ 
  estimation.
- **Converts heterogeneous units** — flux (mmol gDCW⁻¹ h⁻¹) and proteomics (ppm or 
  mol gDCW⁻¹) are reconciled to produce $k_\text{app}$ in s⁻¹.
- **Resolves substrates** — filters cofactors and identifies the acting substrate per 
  reaction.
- **Aggregates over conditions** — computes $k_\text{max}$ (the maximum $k_\text{app}$ 
  across conditions) as a lower-bound estimate of the enzyme's in vivo catalytic 
  capacity, following Davidi et al. 2016.
- **Benchmarks against in vitro data** — links results to EnzyExtractDB (Wei et al. 
  2025) via enzyme sequence and substrate SMILES for *in vivo* / *in vitro* 
  comparison.


## Navigation

* [Getting Started](#getting-started)
* [Basic usage](#basic-usage)
* [Repository Structure](#repository-structure)


## Getting Started

1. Clone the repository:
   ```bash
   git clone https://github.com/lyach/BactoCat.git
   ```

2. Navigate to the project directory:
   ```bash
   cd BactoCat
   ```

3. Create an isolated virtual environment:
   
   **Using `uv` (recommended)**
    ```bash
    # Install uv (if first time)
   curl -LsSf https://astral.sh/uv/install.sh | sh # Linux/macOS
   # or
   pip install uv # Windows

   # Create environment
   uv venv --python 3.10 .venv

   # Activate environment
   source .venv/bin/activate # Linux/macOS
   # or 
   .venv\Scripts\activate # Windows

   # Install dependencies
   uv sync
   ```

   **Using `conda`**
    ```bash
    # Create environment
   conda create -n bactocat python=3.10

   # Activate environment
   conda activate bactocat

   # Install dependencies
   pip install -e .
   ```


## Basic Usage

After set-up, you can run the main pipeline  with:
```bash
python scripts/run_kapp_pipeline.py configs/run_kapp_pipeline/davidi_consensus.yaml
```

You can also run a detailed verbose output with `v`:
```bash
python scripts/run_kapp_pipeline.py -v configs/run_kapp_pipeline/davidi_consensus.yaml
```

A full documentation of the pipeline parameters, inputs and outputs can be found [here](docs/pipeline_specs.md).


## Repository Structure

```

├── configs            <- Configuration files for automated runs
│
├── data
│   ├── external       <- Data from external pipelines
│   ├── interim        <- Intermediate data
│   ├── processed      <- Canonical data ready for modeling
│   └── raw            <- Original data
│
├── docs               <- Extended documentation
│
├── results            <- Store pipeline outputs
│
├── scripts            <- Run pipeline & analyses
│
└── src                           <- Source code for BactoCat
    │
    ├── __init__.py               <- Makes a Python module
    │
    ├── config.py                 <- Stores variables 
    │
    ├── enzyme_classifier.py      <- Functions for enzyme classification using GEM GPR rules
    │
    ├── gene_sequence_maper.py    <- Functions to recover protein sequences via UniProt
    │
    ├── kap_builder.py            <- Functions to build kapp and kmax datasets
    │
    ├── paxdb_mapper.py           <- Functions to map gene IDs to PaxDB proteomics data
    │
    ├── plots.py                  <- Functions to create visualizations
    │
    ├── substrate_mapper.py       <- Functions to recover enzyme substrates from GEM reactions
    │ 
    └── utils.py                  <- Utilities

```

## Project Status

**BactoCat** is a project under active development. The current release implements the pipeline using a published dataset of  *E. coli* growth and proteomics across 31 conditions (Schmidt et al; 2016).