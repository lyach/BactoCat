# k<sub>app</sub> Pipeline

This directory contains an automated pipeline for calculating apparent *in vivo* turnover numbers (k<sub>app</sub>) for *E. coli* enzymes.



## Quick Start

### Basic Usage

```bash
cd BactoCat
python .\scripts\run_kapp_pipeline.py .\configs\scripts\ecoli_homomeric.yaml
```

## Configuration File

The pipeline is controlled by a YAML configuration file (saved under `configs\scripts`). Here's what each parameter does:

### Model Configuration

- **`organism`**: Name of the organism, to use in outputs
- **`model_path`**: Path to your SBML model file (e.g., `iML1515_GEM.xml`)
- **`flux_method`**: Choose between `"FBA"` or `"pFBA"` for flux simulations

### Flux Simulation Parameters

- **`carbon_uptake`**: List of carbon uptake rates in mmol/gDW·h (e.g., `[2, 6, 10]`)
- **`oxygen_uptake`**: List of oxygen uptake rates in mmol/gDW·h (e.g., `[15, 17.5, 20]`)

The pipeline will test all combinations of these rates (e.g., 3 carbon × 3 oxygen = 9 conditions).

### Proteomics Parameters

- **`p_total`**: Total cellular protein content values in g/gDCW (e.g., `[0.32, 0.435, 0.55]`)
- **`paxdb_path`**: Path to PaxDB proteomics data file

### Input Data

If available, one can pass cached data for:

- **`substrate_df`**: Path to substrate information CSV (SMILES and reaction mapping)
- **`sequence_df`**: Path to protein sequence CSV (from UniProt)

### Output

- **`output_dir`** *(optional)*: Directory for results. Defaults to `../scripts/results/`

## Output Files

The pipeline generates:

- **`iml1515_homomeric_kmax_{method}_variability.csv`**: Main results file containing:
  - Maximum k<sub>app</sub> values for each enzyme-substrate pair
  - η statistics (mean, standard deviation, CV)
  - Metadata (genes, reactions, subsystems)


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

### Import Errors

If you see import errors for `enzyme_classifier` or `kapp_builder`, ensure:
1. You're running the script from the `experiments/` directory
2. The `scripts/` directory exists and contains the required modules

### File Not Found Errors

- Check that all paths in `config.yaml` are relative to the `experiments/` directory
- Verify that input files exist at the specified locations

### Optimization

Large GEMs and pFBA most likely require full license solvers (CPLEX or Gurobi). glpk 

## Advanced Usage

### Multiple Configurations

You can create multiple config files for different analyses:

```bash
python run_kapp_pipeline.py ecoli_homomeric.yaml
python run_kapp_pipeline.py pputida_homomeric.yaml
python run_kapp_pipeline.py bsubtilis_homomeric.yaml
```

### Programmatic Usage

You can also import and use the pipeline function directly:

```python
from run_kapp_pipeline import run_kapp_pipeline

kapp_results, kmax_results = run_kapp_pipeline(
    model_path="../",
    flux_method="pFBA",
    carbon_uptake=[2, 6, 10],
    oxygen_uptake=[15, 17.5, 20],
    p_total=[0.32, 0.435, 0.55],
    substrate_df="../",
    sequence_df="../",
    paxdb_path="../"
)
```


