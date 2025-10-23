# BactoCat Dataset

A comprehensive dataset and analysis pipeline for *in vivo*-like kcat (catalytic turnover number) values, integrating multiple databases and metabolic models for enzyme kinetics research.

## Overview

This repository contains curated datasets, processing scripts, and analysis notebooks for studying enzyme kinetics across different organisms.


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

**Updating the environment**

Whenever you modify dependencies (add or update a package), run:
```
uv sync --upgrade
```

or to install a single new package:
```
uv add <package-name>
```


## Solver set-up

The BactoCat pipeline relies on the IBM CPLEX Optimization Studio solver.
To use CPLEX within this environment, install the Python API from your local CPLEX Studio installation (requires a valid academic or commercial license):
```
pip install "C:\Program Files\IBM\ILOG\CPLEX_Studio2211\cplex\python\3.10\x64_win64"
```

  - **NOTE**: CPLEX is not included on `environment.yml`. Each user must install it manually from their own licensed CPLEX Studio installation. Ensure your license environment variables (e.g., ILOG_LICENSE_FILE) are correctly configured.



## Repository Structure

### 📁 `data/`
Main data directory organized by processing stage. Only the `/raw` and `/results` folder are kept out of the `.gitignore`.

- `data/raw/` - Original, unprocessed datasets

- `data/processed/` - Processed and cleaned datasets
  - TO DO: Add subfolders

- `data/final/` - Final curated datasets

- `data/discarded/` - Archived datasets not used in final analysis. Contains previous versions of cleaned and processed files that didn't make it to final analysis.

### 📁 `scripts/`
Data processing and analysis scripts:
- **`enzyme_classifier.py`** - Script for classifying enzymes and their properties
- TO DO: Add the rest

### 📁 `experiments/`
- **`ECOMICS_fluxomics.ipynb`** - Complete processing and analysis of ECOMICS fluxomics data for *E. coli* 
- **`flux_reconciliation.ipynb`** - Pipeline for reconciling ECOMICS fluxes using a 2 level LP approach. Discarded because of poor performance with ECOMICS data.
- **`kapp_pipeline.ipynb`** - Complete processing and analysis of *in vivo*
 $k_{cat}$ for *E. coli* 
- **`load_invitro_datasets.ipynb`** - Loading and analyses of *in vitro* $k_{cat}$ for *E. coli* 
- **`tsne_invivo_invitro.ipynb`** - Analysis of the functional space of *in vivo* and *in vitro* $k_{cat}$ for *E. coli* 

### 📁 `old_experiments/`, 📁 `old_scripts/`
Archived experimental notebooks and scripts

### 📁 `results/`
Generated results from analyses

### 📁 `misc/`
Miscellaneous files


### 📁 `QE/`
Scripts used in Qualifier Examination results