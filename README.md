# CatNivo Dataset

A comprehensive dataset and analysis pipeline for *in vivo*-like kcat (catalytic turnover number) values, integrating multiple databases and metabolic models for enzyme kinetics research.

## Overview

This repository contains curated datasets, processing scripts, and analysis notebooks for studying enzyme kinetics across different organisms.

## Usage

1. **Data exploration**: Start with notebooks in `experiments/` for current analyses
2. **Data processing**: Use scripts in `scripts/` for data manipulation
3. **Results**: Find processed outputs in `data/final`

## Environment set-up

This project uses a reproducible Conda environment specified in `environment.yml`

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
conda env update -f environment.yml --prune
```
- *NOTE:* `prune` removes packages not listed in the YAML. Do not include in the command if you've installed and wish to keep additional packages.


## Repository Structure

### 📁 `data/`
Main data directory organized by processing stage:

#### `data/raw/` - Original, unprocessed datasets

#### `data/processed/` - Processed and cleaned datasets
- **`ECOMICS/`** - Processed ECOMICS data
  - `fluxomics_ecomics_mapped.csv` - Mapped fluxomics data
  - `eco_out_misses_curated.csv` - Curated missing data
  - `iml1515_rxns.csv` - iML1515 model reactions
- **`UniProt/`** - Protein mapping files from UniProt queries

#### `data/final/` - Final curated datasets
- Currently empty - reserved for final processed datasets

#### `data/discarded/` - Archived datasets not used in final analysis
- Contains various versions of cleaned and processed files that didn't make it to final analysis

### 📁 `scripts/`
Data processing and analysis scripts:
- **`enzyme_classifier.py`** - Script for classifying enzymes and their properties

### 📁 `experiments/`
- **`ECOMICS_fluxomics.ipynb`** - Analysis of ECOMICS fluxomics data

### 📁 `old_experiments/`
Archived experimental notebooks

### 📁 `old_scripts/`
Archived processing scripts

### 📁 `results/`
Generated analysis results and outputs

### 📁 `misc/`
Miscellaneous files


