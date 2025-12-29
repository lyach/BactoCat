# BactoCat 

A computational framework designed to bridge the gap between *in vitro* and *in vivo* enzyme kinetics. 

By integrating large-scale phenomics, proteomics, and genome-scale metabolic models (GEMs), BactoCat derives apparent catalytic rates ($k_{app}$) that reflect the true catalytic potential of enzymes within bacterial cells.

## Table of Contents
* [Getting Started](#getting-started)
* [Basic usage](#basic-usage)
* [Environment](#environment-set-up)
* [Solver](#solver-set-up)
* [Repository Structure](#repository-structure)

---

## Getting Started

Clone the repository  :
```bash
git clone https://github.com/lyach/BactoCat.git
```

Once cloned, you can install BactoCat as an editable package:

```bash
cd BactoCat
pip install -e .
```

This installs:
- All required dependencies
- The `run-kapp-pipeline` CLI command
- BactoCat modules for import


## Basic Usage

After installation, you can run the pipeline using the CLI command:

```bash
run-kapp-pipeline configs/run_kapp_pipeline/ecoli_homomeric.yaml
```

Or run the module using Python or uv:
```bash
# python
python -m scripts.run_kapp_pipeline configs/run_kapp_pipeline/ecoli_homomeric.yaml

# uv 
uv run -m scripts.run_kapp_pipeline configs/run_kapp_pipeline/ecoli_homomeric.yaml
```

You can also run a detailed verbose output with:
```bash
run-kapp-pipeline -v configs/run_kapp_pipeline/ecoli_homomeric.yaml
```

A detailed documentation of the pipeline parameters can be found [here](docs/kapp_pipeline.md).

## Environment set-up

### Using `conda`

A reproducible Conda environment can be built using the `environment.yml`

**1. Create the environment (first time)**
```
conda env create -f environment.yml
```

**2. Activate the environment**
```
conda activate bactocat_env
```

**3. Update the environment**
```
conda env update -f environment.yml 
```
- *NOTE:* Include `--prune` at the end of this command to remove packages not listed in the YAML. If you've installed and wish to keep additional packages, run as it is.

### Using `uv`

The following steps will set up a reproducible Python environment using `uv`. This approach avoids dependency conflicts and works seamlessly on HPC clusters.

**1. Install uv**

This will install a single static binary in your `$HOME/bin`:
```
curl -LsSf https://astral.sh/uv/install.sh | sh
```
**2. Navigate to the repository**
```
cd ~/BactoCat
```

**3. Create a virtual environment**

Creates an isolated environment inside the repository:
```
uv venv .venv
```

**4. Activate venv**
```
source .venv/bin/activate
```

**5. Install dependencies**
Synchronize and install all dependencies listed in `pyproject.toml`:
```
uv sync
```
   This automatically resolves and locks dependencies for reproducibility (`uv.lock` file).

**6. Updating the environment**

Whenever you modify dependencies (add or update a package), run:
```
uv sync --upgrade
```

To install a single new package:
```
uv add <package-name>
```

To remove a single new package:
```
uv remove <package-name>
``` 


## Solver set-up

The BactoCat framework relies on external solvers to perform optimization tasks like parsimonious Flux Balance Analysis (pFBA). To run the **complete pipeline**, you must have at least one of the following solvers installed on your system:
- [IBM CPLEX](https://www.ibm.com/products/ilog-cplex-optimization-studio) (requires license) 
- [Gurobi](https://www.gurobi.com/lp/all/request-an-evaluation/?utm_source=google&utm_medium=cpc&utm_campaign=M3+DG+Search+NA+Brand&gad_source=1&gad_campaignid=193283256&gbraid=0AAAAADimQ3hP-_zFR-f3IsRZ0E_dQ7VmZ&gclid=CjwKCAiAmKnKBhBrEiwAaqAnZy5G51otRX82XLRp-6SEAtRUtVB540b6DdUNlO04R-mmwWGvxJxVHxoCBKgQAvD_BwE) (requires license)

**Note**: The `cplex` and `gurobipy` packages are not included in the dependencies. You must manually install your licensed executable.

**Alternative**: If you do not have a commercial solver license, you can run **Light-BactoCat**. This configuration uses the `glpk` solver for optimization. Note that Flux Variability Analysis (FVA) is disabled in this mode.

## Repository Structure

### 📁 `data/`
Main data directory organized by processing stage. Only the `/raw` and `/results` folder are kept out of the `.gitignore`.

- `data/raw/` - Original, unprocessed datasets and genome-scale models

- `data/processed/` - Processed datasets, output from notebooks and scripts
- `data/final/` - Final datasets


### 📁 `src/`
Main source code for the BactoCat pipelines:
- **`enzyme_classifier.py`** - Script for classifying enzymes and their properties
- TO DO: Add the rest

### 📁 `notebooks/`
- **`ECOMICS_fluxomics.ipynb`** - Complete processing and analysis of ECOMICS fluxomics data for *E. coli* 
- **`flux_reconciliation.ipynb`** - Pipeline for reconciling ECOMICS fluxes using a 2 level LP approach. Discarded because of poor performance with ECOMICS data.
- **`kapp_pipeline.ipynb`** - Complete processing and analysis of *in vivo*
 $k_{cat}$ for *E. coli* 
- **`load_invitro_datasets.ipynb`** - Loading and analyses of *in vitro* $k_{cat}$ for *E. coli* 
- **`tsne_invivo_invitro.ipynb`** - Analysis of the functional space of *in vivo* and *in vitro* $k_{cat}$ for *E. coli* 

### 📁 `scripts/`
Scripts to run BactoCat pipelines automated with config files 

### 📁 `misc/`
Miscellaneous and discarded files