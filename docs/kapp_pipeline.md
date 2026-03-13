# k<sub>app</sub> Pipeline


## Structure

- The main BactoCat pipeline is executed through the [`run_kapp_pipeline`](../scripts/run_kapp_pipeline.py) script.
- It requires a `YAML` configuration file specifying model parameters, flux simulation conditions, proteomics data, and filtering thresholds. 
- It will processes genome-scale metabolic models through multiple steps, producing curated datasets of apparent turnover numbers (k<sub>app</sub>).


## Configuration File

The pipeline is controlled by a `YAML` configuration file (saved under `configs/run_kapp_pipeline/`). All paths are relative to the **project root**.

### Model Configuration

- **`organism`**: Organism name for outputs
- **`model_path`**: Path to SBML model file 
- **`solver`**: Optimization solver - `"cplex"`, `"gurobi"`, or `"glpk"`
- **`flux_method`**: Flux analysis method - `"FBA"` or `"pFBA"`

### Flux Simulations

The pipeline can be run with ways of computng flux simulations:

#### a) Using experimental conditions:

- **`medium_df`**: csv file specifying the uptake reactions (columns) across different conditions (rows). Values are in  in mmol/gDW·h 

An example file can be found [here](../configs/run_kapp_pipeline/ecoli_medium_aidaying.yaml).

#### b) Using different combinations of uptake rates:

- **`carbon_uptake`**: Carbon uptake rates in mmol/gDW·h 
- **`oxygen_uptake`**: Oxygen uptake rates in mmol/gDW·h 
- **`carbon_exchange_rxn`**: Carbon exchange reaction ID 
- **`oxygen_exchange_rxn`**: Oxygen exchange reaction ID
- **`mu_fraction`**: Growth rate fraction for FVA 

The pipeline tests all combinations of carbon and oxygen rates. An example file can be found [here](../configs/run_kapp_pipeline/ecoli_homomeric.yaml).

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

### Output folder
- **`folder_id`**: name of the folder to store the run results.


## Output Files

Results will be saved under `results/run_kapp_pipeline/{folder_id}/`:

```
results/run_kapp_pipeline/folder_id/
├── data/                                 # Intermediate files (sequences and substrates)
└── results/                              # Final outputs
    ├── kmax.csv                          # Main results
    ├── log_ecoli_{folder_id}.log         # Detailed execution log
    ├── FBA_fluxomics.parquet             # All fluxomes across conditions
    └── FVA_violations.csv                # Fluxes violating FVA
```

### Main Results File

`kmax.csv` contains:
- Maximum k<sub>app</sub> values for each enzyme-substrate pair
- η (saturation) statistics: mean, standard deviation, coefficient of variation
- Metadata: genes, reactions, subsystems, sequences, SMILES


## Pipeline Steps

The automated pipeline performs the following steps:

1. **Load Model**: Read SBML metabolic model
2. **Extract GPR Rules**: Classify enzymes (homomeric, complex, isoenzyme)
3. **Flux Simulations**: Run FBA/pFBA for all conditions
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

### Solvers

The BactoCat framework relies on external solvers to perform optimization tasks. To run certain components like FVA, you must have at least one of the following solvers installed on your system:
- [IBM CPLEX](https://www.ibm.com/products/ilog-cplex-optimization-studio) (requires license) 
- [Gurobi](https://www.gurobi.com/lp/all/request-an-evaluation/?utm_source=google&utm_medium=cpc&utm_campaign=M3+DG+Search+NA+Brand&gad_source=1&gad_campaignid=193283256&gbraid=0AAAAADimQ3hP-_zFR-f3IsRZ0E_dQ7VmZ&gclid=CjwKCAiAmKnKBhBrEiwAaqAnZy5G51otRX82XLRp-6SEAtRUtVB540b6DdUNlO04R-mmwWGvxJxVHxoCBKgQAvD_BwE) (requires license)

**Notes**: 
- The `cplex` and `gurobipy` packages are not included in the dependencies. You must manually install your licensed executable.
- If you do not have a commercial solver license, you can still run the pipeline, using the default `glpk` solver in the COBRApy package. Note that Flux Variability Analysis (FVA) cannot be run.

