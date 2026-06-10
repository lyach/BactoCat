# BactoCat Pipeline


## Structure

- The main BactoCat pipeline is executed through the [`run_kapp_pipeline`](../scripts/run_kapp_pipeline.py) script.
- It requires a `YAML` configuration file specifying model parameters.
- It processes genome-scale metabolic models through multiple steps, producing curated datasets of maximal (k<sub>max</sub>) and apparent (k<sub>app</sub>) *in vivo* turnover numbers.

### Usage

```bash
run-kapp-pipeline configs/run_kapp_pipeline/ecoli_homomeric.yaml
python -m scripts.run_kapp_pipeline configs/run_kapp_pipeline/ecoli_homomeric.yaml
```

Add `-v` / `--verbose` for DEBUG-level console output.


## Configuration File

The pipeline is controlled by a `YAML` configuration file (saved under `configs/run_kapp_pipeline/`). All paths are relative to the **project root**.

### Model Configuration

| Parameter | Type | Required | Default | Description |
|---|---|---|---|---|
| `organism` | str | ✓ | — | Organism name (e.g. `"ecoli"`) |
| `model_path` | path | ✓ | — | Path to SBML model file (`.xml`) |
| `folder_id` | str | | `"run"` | ID appended to organism name for the output directory |
| `solver` | str | | `"cplex"` | Optimization solver: `"cplex"`, `"gurobi"`, or `"glpk"` |

### Flux Simulations

| Parameter | Type | Required | Default | Description |
|---|---|---|---|---|
| `flux_method` | str | | `"pFBA"` | Flux analysis method: `"FBA"` or `"pFBA"` |
| `medium_df` | path | ✓ | — | CSV specifying uptake reactions (columns) across conditions (rows), in mmol/gDW·h |
| `include_growth` | bool | | `false` | Fix growth reaction bounds to the measured value in `medium_df` |
| `free_metabolites` | list[str] | | `null` | Exchange reaction IDs to leave unconstrained during simulations |
| `carbon_uptake` | list[float] | | `null` | Carbon uptake rates to test (mmol/gDW/h); alternative to `medium_df` grid |
| `oxygen_uptake` | list[float] | | `null` | Oxygen uptake rates to test (mmol/gDW/h); alternative to `medium_df` grid |
| `carbon_exchange_rxn` | str | | `"EX_glc__D_e"` | Reaction ID for carbon exchange |
| `oxygen_exchange_rxn` | str | | `"EX_o2_e"` | Reaction ID for oxygen exchange |
| `run_fva` | bool | | `true` | Run FVA to filter thermodynamically infeasible fluxes |
| `mu_fraction` | float | | `0.9` | Fraction of optimal growth rate used as FVA lower bound |

### Proteomics Parameters

Two modes are supported, controlled by `specific_proteome`:

**Consensus proteome** (`specific_proteome: false`, default):

| Parameter | Type | Required | Default | Description |
|---|---|---|---|---|
| `specific_proteome` | bool | | `false` | Use condition-matched proteomics instead of PaxDB |
| `paxdb_path` | path | ✓ | — | Path to PaxDB proteomics data file |
| `p_total` | float | ✓ | — | Total protein fraction (g protein / g DCW) |

**Condition-matched proteome** (`specific_proteome: true`):

| Parameter | Type | Required | Default | Description |
|---|---|---|---|---|
| `specific_proteome` | bool | ✓ | — | Set to `true` to use condition-matched data |
| `proteomics_path` | path | ✓ | — | Path to condition-matched proteomics CSV |

### Input Data (Optional)

Pre-computed data can be provided to skip auto-generation:

| Parameter | Type | Default | Description |
|---|---|---|---|
| `substrate_df` | path | `null` | Path to substrate CSV (SMILES and reaction mapping); auto-generated if absent |
| `sequence_df` | path | `null` | Path to protein sequence CSV (from UniProt); auto-generated if absent |

### Filtering Thresholds

| Parameter | Type | Default | Description |
|---|---|---|---|
| `upper_threshold` | float | `1.0e6` | Maximum k<sub>app</sub> (s⁻¹) |
| `lower_threshold` | float | `1.0e-5` | Minimum k<sub>app</sub> (s⁻¹) |

### Output Options

| Parameter | Type | Default | Description |
|---|---|---|---|
| `calculate_eta` | bool | `true` | Compute η = k<sub>app</sub> / k<sub>max</sub> and variance metrics across conditions |
| `save_kapp` | bool | `false` | Save the full wide-format k<sub>app</sub> dataframe as `kapp.csv` |


## Output Files

The run name is constructed as `{organism}_{folder_id}`. Results are saved under `results/run_kapp_pipeline/{run_name}/`:

```
results/run_kapp_pipeline/{organism}_{folder_id}/
├── data/
│   ├── {organism}_uniprot_seqs.csv          # Protein sequences (auto-generated only)
│   └── {organism}_substrate_df.csv          # Substrate SMILES (auto-generated only)
└── results/
    ├── kmax.csv                              # kmax values
    ├── kapp.csv                              # Wide-format kapp
    └── log_{run_name}.log                    # Detailed execution log
```

### `kmax.csv`

Contains, for each enzyme:
- k<sub>max</sub> value (maximum k<sub>app</sub> across conditions)
- Metadata: genes, reactions, subsystems, protein sequences, substrate SMILES
- η statistics (mean, variance, coefficient of variation) if `calculate_eta=true`


## Pipeline Steps

```
Model setup (solver + SBML load)
│
├── STEP 1  — Extract GPR rules: classify enzymes (homomeric, complex, isoenzyme)
├── STEP 2  — Run FBA/pFBA flux simulations across all conditions
├── STEP 3  — [Optional] Run FVA and filter infeasible fluxes
├── STEP 4  — Load protein sequences (from file or UniProt API)
├── STEP 5  — Load substrate information (from file or model extraction)
├── STEP 6  — Merge enzyme, flux, substrate, and sequence data
├── STEP 7  — Map proteomics (PaxDB consensus or condition-matched)
├── STEP 8  — Calculate k<sub>app</sub> for homomeric enzymes
├── STEP 9  — Filter k<sub>app</sub> values outside physical thresholds
├── STEP 10 — Determine k<sub>max</sub> (maximum k<sub>app</sub> per enzyme across conditions)
├── STEP 11 — [Optional] Calculate η = k<sub>app</sub> / k<sub>max</sub> per condition
└── STEP 12 — Save results (kmax.csv, and kapp.csv if save_kapp=true)
```


## Troubleshooting

### Solvers

The BactoCat framework relies on external solvers to perform optimization tasks. To run certain components like FVA, you must have at least one of the following solvers installed:

- [IBM CPLEX](https://www.ibm.com/products/ilog-cplex-optimization-studio) (requires license)
- [Gurobi](https://www.gurobi.com) (requires license)

**Notes**:
- The `cplex` and `gurobipy` packages are not included in the dependencies. You must manually install your licensed executable.
- If you do not have a commercial solver license, you can still run the pipeline using the default `glpk` solver bundled with COBRApy. Note that Flux Variability Analysis (`run_fva: true`) requires CPLEX or Gurobi.
