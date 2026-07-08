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

### Enzyme Classes

| Parameter | Type | Default | Description |
|---|---|---|---|
| `enzyme_classes` | list[str] | `["homomeric"]` | Enzyme classes to compute k<sub>app</sub>/k<sub>max</sub> for. Any subset of `homomeric`, `isoenzyme`, `complex`, `mixed`. |

- `homomeric`: single-gene enzymes; k<sub>app</sub> = (flux/3600) / abundance (s⁻¹).
- `isoenzyme` / `complex` / `mixed`: computed with the Davidi mass-fraction / specific-activity method (assuming active sites = 1 and subunit stoichiometry = 1). `mixed` covers isozymes that are themselves complexes. See [`computation.md`](computation.md) for the equations, assumptions, and missing-data policy.

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

Contains, for each enzyme (homomeric protein or integrated isoenzyme/complex):
- `kcat_app_max` — k<sub>max</sub>, the maximum turnover k<sub>app</sub> (s⁻¹) across conditions, with `condition_max`
- `SA_app_max` — maximum specific activity (µmol·mg⁻¹·min⁻¹) across conditions
- `enzyme_class` (`homomeric`, `isoenzyme`, `complex`, `mixed`), `enzyme_key` (grouping identity)
- `mass_fraction` (mg/gDCW) and `MW_total` (g/mol) used in the denominator
- `n_components` (genes in the GPR) and `n_contributing` (genes summed after the missing-data policy)
- Metadata: `genes`, `sequence`/`sequences`, `rxn`, `subsystem`, substrate `SMILES`, `Direction`
- η statistics (mean, variance, coefficient of variation) if `calculate_eta=true`

For homomeric enzymes `SA_app_max` corresponds to the same condition as k<sub>max</sub>; for grouped enzymes whose contributing subunit set varies across conditions the two maxima may come from different conditions.


## Pipeline Steps

```
Model setup (solver + SBML load)
│
├── STEP 1  — Extract GPR rules: DNF-parse and classify enzymes (homomeric, complex, isoenzyme, mixed)
├── STEP 2  — Run FBA/pFBA flux simulations across all conditions
├── STEP 3  — [Optional] Run FVA and filter infeasible fluxes
├── STEP 4  — Load protein sequences (from file or UniProt API)
├── STEP 5  — Load substrate information (from file or model extraction)
├── STEP 6  — Merge enzyme, flux, substrate, sequence data; compute molecular weights
├── STEP 7  — Map proteomics (PaxDB consensus or condition-matched)
├── STEP 8  — Calculate k<sub>app</sub> for the requested `enzyme_classes` (homomeric turnover; isoenzyme/complex/mixed via mass-fraction SA + s⁻¹)
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
