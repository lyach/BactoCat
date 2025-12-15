# k<sub>app</sub> Pipeline

This directory contains an automated pipeline for calculating apparent *in vivo* turnover numbers (k<sub>app</sub>) from genome-scale metabolic models.

## Installation

First, install BactoCat as an editable package from the project root:

```bash
cd BactoCat
pip install -e .
```

This installs:
- All required dependencies
- The `run-kapp-pipeline` CLI command
- BactoCat modules for import

## Quick Start

### Basic Usage

After installation, run the pipeline using the CLI command:

```bash
# From anywhere (after installation)
run-kapp-pipeline configs/run_kapp_pipeline/ecoli_homomeric.yaml

# Or using Python module syntax
python -m scripts.run_kapp_pipeline configs/run_kapp_pipeline/ecoli_homomeric.yaml
```

### Verbose Mode

For detailed debug output including function calls and line numbers:

```bash
run-kapp-pipeline -v configs/run_kapp_pipeline/ecoli_homomeric.yaml
```

## Configuration File

The pipeline is controlled by a YAML configuration file (saved under `configs/run_kapp_pipeline/`). All paths are relative to the **project root**.

### Model Configuration

- **`organism`**: Organism name for outputs
- **`model_path`**: Path to SBML model file 
- **`solver`**: Optimization solver - `"cplex"`, `"gurobi"`, or `"glpk"`
- **`flux_method`**: Flux analysis method - `"FBA"` or `"pFBA"`

### Flux Simulation Parameters

- **`carbon_uptake`**: Carbon uptake rates in mmol/gDW·h 
- **`oxygen_uptake`**: Oxygen uptake rates in mmol/gDW·h 
- **`carbon_exchange_rxn`**: Carbon exchange reaction ID 
- **`oxygen_exchange_rxn`**: Oxygen exchange reaction ID
- **`mu_fraction`**: Growth rate fraction for FVA 

The pipeline tests all combinations of carbon and oxygen rates (e.g., 3 × 3 = 9 conditions).

### Proteomics Parameters

- **`p_total`**: Total protein fractions in g/gDCW 
- **`paxdb_path`**: Path to PaxDB proteomics data file

### Input Data (Optional)

Pre-computed data can be provided to speed up the pipeline:

- **`substrate_df`**: Path to substrate CSV (SMILES and reaction mapping)
- **`sequence_df`**: Path to protein sequence CSV (from UniProt)

If not provided, these will be auto-generated.

### Filtering Thresholds

- **`upper_threshold`**: Maximum k<sub>app</sub> in s⁻¹ 
- **`lower_threshold`**: Minimum k<sub>app</sub> in s⁻¹

### Example Configuration

```yaml
# configs/run_kapp_pipeline/ecoli_homomeric.yaml
organism: "ecoli"
model_path: "data/raw/gems/iml1515.xml"
solver: "cplex"
flux_method: "pFBA"

carbon_uptake: [2, 6, 10]
oxygen_uptake: [15, 17.5, 20]
carbon_exchange_rxn: "EX_glc__D_e"
oxygen_exchange_rxn: "EX_o2_e"
mu_fraction: 0.9

p_total: [0.32, 0.435, 0.55]
paxdb_path: "data/raw/paxdb/ecoli_whole_organism_integrated_511145.txt"

upper_threshold: 1.0e6
lower_threshold: 1.0e-5

sequence_df: "data/processed/uniprot/ecoli_uniprot_seqs.csv"
substrate_df: "data/processed/substrates/iml1515_substrates.csv"

```

## Output Files

Results are saved to `results/run_kapp_pipeline/{organism}_{date}_{id}/`:

```
results/run_kapp_pipeline/ecoli_20251215_1634/
├── data/                                        # Intermediate files (sequences and substrates)
└── results/                                        # Final outputs
    ├── kmax_ecoli_20251215_1634.csv                # Main results
    ├── log_ecoli_20251215_1634.log                 # Detailed execution log
    ├── FVA_bounds_ecoli_20251215_1634.csv          # FVA min/max bounds
    ├── fluxomics_filtered_ecoli_20251215_1634.csv  # FVA-filtered fluxes
    └── FVA_violations_ecoli_20251215_1634.csv      # Fluxes violating FVA
```

### Main Results File

`kmax_{organism}_{date}_{id}.csv` contains:
- Maximum k<sub>app</sub> values for each enzyme-substrate pair
- η (saturation) statistics: mean, standard deviation, coefficient of variation
- Metadata: genes, reactions, subsystems, sequences, SMILES


## Pipeline Steps

The automated pipeline performs the following steps:

1. **Load Model**: Read SBML metabolic model
2. **Extract GPR Rules**: Classify enzymes (homomeric, complex, isoenzyme)
3. **Flux Simulations**: Run FBA/pFBA for all condition combinations
4. **Load Sequences**: Import protein sequences from UniProt
5. **Load Substrates**: Import substrate SMILES information
6. **Merge Data**: Combine enzyme, flux, substrate, and sequence data
7. **Map Proteomics**: Link protein abundance data from PaxDB
8. **Calculate k<sub>app</sub>**: Compute apparent turnover numbers
9. **Filter Values**: Remove unrealistic k<sub>app</sub> values
10. **Determine k<sub>max</sub>**: Find maximum k<sub>app</sub> across conditions
11. **Calculate η**: Compute enzyme saturation metrics 
12. **Save Results**: Export final dataframe



## Troubleshooting

### Command Not Found

If `run-kapp-pipeline` is not recognized:
1. Ensure you've run `pip install -e .` from the project root
2. Check that your virtual environment is activated
3. Use `python -m scripts.run_kapp_pipeline` as an alternative

### Import Errors

If you see `ModuleNotFoundError`:
```bash
# Reinstall the package
pip install -e .
```

### Configuration Validation Errors

If you see "Model file not found" or similar errors:
- **All paths in YAML must be relative to the project root**, not the YAML file location

### Solver Issues

- Large GEMs and pFBA require commercial solvers (CPLEX or Gurobi). Ensure your solver is properly licensed and accessible



## Advanced Usage

### Multiple Configurations

Create separate config files for different organisms or conditions:

```bash
run-kapp-pipeline configs/run_kapp_pipeline/ecoli_homomeric.yaml
run-kapp-pipeline configs/run_kapp_pipeline/pputida_homomeric.yaml
run-kapp-pipeline configs/run_kapp_pipeline/yeast_homomeric.yaml
```

### Programmatic Usage

Import and use the pipeline function directly in Python scripts:

```python
from pathlib import Path
from src.config import PipelineConfig
from scripts.run_kapp_pipeline import run_kapp_pipeline

# Create configuration programmatically
config = PipelineConfig(
    organism="ecoli",
    model_path=Path("data/raw/gems/iml1515.xml"),
    solver="cplex",
    flux_method="pFBA",
    carbon_uptake=[2, 6, 10],
    oxygen_uptake=[15, 17.5, 20],
    carbon_exchange_rxn="EX_glc__D_e",
    oxygen_exchange_rxn="EX_o2_e",
    p_total=[0.32, 0.435, 0.55],
    paxdb_path=Path("data/raw/paxdb/ecoli_whole_organism_integrated_511145.txt")
)

# Set up output directories
output_dir = Path("results/run_kapp_pipeline/my_run/results")
data_dir = Path("results/run_kapp_pipeline/my_run/data")
output_dir.mkdir(parents=True, exist_ok=True)
data_dir.mkdir(parents=True, exist_ok=True)

# Run pipeline
kapp_results, kmax_results = run_kapp_pipeline(
    config=config,
    output_dir=output_dir,
    data_dir=data_dir,
    run_name="my_run"
)
```

