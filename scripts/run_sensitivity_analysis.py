"""
BactoCat Sensitivity Analysis Runner.

Runs the kapp pipeline at multiple total protein contents (p_total)
to assess sensitivity of kmax estimates to the assumed protein fraction.

Steps 1-6 of the kapp pipeline are executed once. 
Steps 7-12 (proteomics mapping, kapp/kmax calculation) are repeated for each p_total.
The resulting kmax values are collated into a single kmax_sensitivity.csv.

Usage:
    python scripts/run_sensitivity_analysis.py configs/run_sensitivity_analysis/your_config.yaml
"""

import argparse
import sys
from datetime import datetime
from pathlib import Path

import cobra
import numpy as np
import pandas as pd
import yaml
from loguru import logger

from src.config import PipelineConfig, load_config, PROJ_ROOT, ensure_dir_exists
from src.enzyme_classifier import create_gpr_dataframe, analyze_model_gprs
from src.gene_sequence_mapper import map_organism_to_uniprot
from src.substrate_mapper import get_substrate_df
from src.proteomics_mapper import paxdb_protein_mapping, specific_proteome_mapping
from src.kapp_builder import (
    create_fluxomics_dataframe,
    create_enzyme_info_dataframe,
    create_FVA_dataframe,
    FVA_integration,
    calculate_kapp_homomeric,
    evaluate_kapp_homomeric,
    get_kmax_homomeric,
    get_eta,
    get_kapp_dataframe,
)


def load_sensitivity_config(yaml_path: str | Path) -> dict:
    """
    Load a sensitivity analysis YAML configuration.

    Parameters
    ----------
    yaml_path : str or Path
        Path to the YAML file containing ``kapp_config``,
        ``p_total_min``, ``p_total_max``, ``p_total_step``,
        and an optional ``folder_id``.

    Returns
    -------
    dict
        Raw configuration dictionary.
    """
    yaml_path = Path(yaml_path)
    if not yaml_path.exists():
        raise FileNotFoundError(f"Configuration file not found: {yaml_path}")
    with open(yaml_path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def setup_logging(log_file: Path, run_name: str, verbose: bool = False) -> None:
    """Configure loguru for console and file output."""
    logger.remove()
    if verbose:
        logger.add(
            sys.stderr,
            format="<green>{time:HH:mm:ss}</green> | <level>{level: <8}</level> | <cyan>{function}</cyan>:<cyan>{line}</cyan> | <level>{message}</level>",
            level="DEBUG",
            colorize=True,
        )
        logger.add(
            log_file,
            format="{time:YYYY-MM-DD HH:mm:ss} | {level: <8} | {function}:{line} | {message}",
            level="DEBUG",
        )
    else:
        logger.add(
            sys.stderr,
            format="<green>{time:HH:mm:ss}</green> | <level>{level: <8}</level> | <level>{message}</level>",
            level="INFO",
            colorize=True,
        )
        logger.add(
            log_file,
            format="{time:YYYY-MM-DD HH:mm:ss} | {level: <8} | {message}",
            level="DEBUG",
            rotation="10 MB",
            retention="30 days",
        )
    logger.info(f"Logging initialized for run: {run_name}")


def run_sensitivity_analysis(
    config: PipelineConfig,
    p_total_values: list[float],
    output_dir: Path,
    data_dir: Path,
    run_name: str,
) -> pd.DataFrame:
    """
    Run the kapp pipeline at multiple p_total values and collate kmax.

    Parameters
    ----------
    config : PipelineConfig
        Validated base pipeline configuration.
    p_total_values : list[float]
        Total protein fractions to sweep (g protein / g DCW).
    output_dir : Path
        Directory for output files.
    data_dir : Path
        Directory for intermediate data files.
    run_name : str
        Name for this run.

    Returns
    -------
    pd.DataFrame
        Merged kmax DataFrame with one ``kcat_app_max_{p}`` column per
        p_total value.
    """
    # ---- Solver ----
    try:
        cobra.Configuration().solver = config.solver
        logger.info(f"Solver set to {config.solver}")
    except Exception as e:
        raise ValueError(f"Error setting solver: {e}")

    # ---- Load model ----
    try:
        model = cobra.io.read_sbml_model(str(config.model_path))
        model = model.copy()
        logger.info(f"Model loaded successfully from {config.model_path.name}")
    except Exception as e:
        raise ValueError(f"Model not found at {config.model_path}. Error: {e}")

    # ==== STEP 1: Create GEM enzymes dataframe ====
    logger.info("=" * 50)
    logger.info("STEP 1: Create GEM enzymes dataframe")
    df_enzymes = create_gpr_dataframe(model)

    stats = analyze_model_gprs(model)
    logger.info(f"Total reactions: {stats['total_reactions']}")
    logger.info(f"Reactions with GPR: {stats['reactions_with_gpr']}")
    logger.info(f"Total genes: {stats['total_genes']}")

    # ==== STEP 2: Run fluxomics simulations ====
    logger.info("=" * 50)
    logger.info(f"STEP 2: Run {config.flux_method} fluxomics simulations")

    medium_df_loaded = None
    if config.medium_df is not None and config.medium_df.exists():
        try:
            medium_df_loaded = pd.read_csv(config.medium_df)
            logger.info(f"Medium dataframe loaded: {len(medium_df_loaded)} conditions from {config.medium_df.name}")
        except Exception as e:
            raise ValueError(f"Error loading medium_df from {config.medium_df}: {e}")

    fluxomics_df = create_fluxomics_dataframe(
        flux_method=config.flux_method,
        GEM=model,
        medium_df=medium_df_loaded,
        include_growth=config.include_growth,
        free_metabolites=config.free_metabolites,
    )

    # ==== STEP 3: Run flux variability analysis ====
    logger.info("=" * 50)
    if not config.run_fva:
        logger.info("STEP 3: Skipping flux variability analysis (run_fva=false)")
        filtered_fluxomics_df = fluxomics_df.copy()
    else:
        logger.info("STEP 3: Run flux variability analysis")
        try:
            fva_df = create_FVA_dataframe(
                GEM_path=str(config.model_path),
                medium_df=medium_df_loaded,
                mu_fraction=config.mu_fraction,
                solver=config.solver,
                include_growth=config.include_growth,
                free_metabolites=config.free_metabolites,
            )
            filtered_fluxomics_df, violations_df = FVA_integration(
                fluxomics_df, fva_df, filter=True
            )
            fva_df.to_csv(output_dir / f"FVA_bounds_{run_name}.csv", index=False)
            filtered_fluxomics_df.to_csv(
                output_dir / f"fluxomics_filtered_{run_name}.csv", index=False
            )
            violations_df.to_csv(
                output_dir / f"FVA_violations_{run_name}.csv", index=False
            )
            logger.info(f"FVA integration complete. Filtered fluxomics: {filtered_fluxomics_df.shape[0]} rows")
        except Exception as e:
            raise RuntimeError(f"Error during FVA integration: {e}")

    fluxomics_df = filtered_fluxomics_df.copy()

    # ==== STEP 4: Load sequence information ====
    logger.info("=" * 50)
    logger.info("STEP 4: Load sequence information")

    if config.sequence_df and config.sequence_df.exists():
        try:
            sequence_df_loaded = pd.read_csv(config.sequence_df)
            logger.info(f"Sequence dataframe loaded: {len(sequence_df_loaded)} rows")
        except Exception as e:
            raise ValueError(f"Sequence dataframe not found at {config.sequence_df}. Error: {e}")
    else:
        logger.info("No sequence dataframe provided, retrieving sequences from UniProt")
        sequence_df_loaded = map_organism_to_uniprot(config.organism)
        seq_output = data_dir / f"{config.organism}_uniprot_seqs.csv"
        sequence_df_loaded.to_csv(seq_output, index=False)
        logger.info(f"Sequence dataframe created: {len(sequence_df_loaded)} rows")

    # ==== STEP 5: Load substrate information ====
    logger.info("=" * 50)
    logger.info("STEP 5: Load substrate information")

    if config.substrate_df and config.substrate_df.exists():
        try:
            substrate_df_loaded = pd.read_csv(config.substrate_df)
            logger.info(f"Substrate dataframe loaded: {len(substrate_df_loaded)} rows")
        except Exception as e:
            raise ValueError(f"Substrate dataframe not found at {config.substrate_df}. Error: {e}")
    else:
        logger.info("No substrate dataframe provided, generating from model")
        substrate_df_loaded = get_substrate_df(model)
        sub_output = data_dir / f"{config.organism}_substrate_df.csv"
        substrate_df_loaded.to_csv(sub_output, index=False)
        logger.info(f"Substrate dataframe created: {len(substrate_df_loaded)} rows")

    # ==== STEP 6: Create enzyme information dataframe ====
    logger.info("=" * 50)
    logger.info("STEP 6: Create enzyme information dataframe")
    enzymes_info_dfs = create_enzyme_info_dataframe(
        df_enzymes, fluxomics_df, substrate_df_loaded, sequence_df_loaded,
        run_fva=config.run_fva,
    )

    # ==== STEPS 7-12: Sweep over p_total values ====
    id_cols = ["sequence", "SMILES", "gene", "rxn", "subsystem"]
    merged: pd.DataFrame | None = None

    for p_total in p_total_values:
        p_label = f"{p_total:.2f}"
        logger.info("=" * 60)
        logger.info(f"SENSITIVITY SWEEP — p_total = {p_label}")
        logger.info("=" * 60)

        # STEP 7: Map proteomics
        logger.info("STEP 7: Map proteomics information")
        if config.specific_proteome:
            enzyme_protein_info_dfs = specific_proteome_mapping(
                enzymes_info_dfs, str(config.proteomics_path)
            )
        else:
            enzyme_protein_info_dfs = paxdb_protein_mapping(
                enzymes_info_dfs, str(config.paxdb_path), p_total=p_total
            )

        # STEP 8: Calculate kapp
        logger.info("STEP 8: Calculate kapp for homomeric enzymes")
        kapp_dfs = calculate_kapp_homomeric(enzyme_protein_info_dfs)

        # STEP 9: Filter
        logger.info("STEP 9: Filter values above physical threshold")
        kapp_dfs_filtered = evaluate_kapp_homomeric(kapp_dfs)

        if config.save_kapp:
            kapp_df = get_kapp_dataframe(kapp_dfs_filtered)
            kapp_df.to_csv(output_dir / f"kapp_p{p_label}.csv", index=False)

        # STEP 10: kmax
        logger.info("STEP 10: Get kmax for homomeric enzymes")
        kmax_df = get_kmax_homomeric(kapp_dfs_filtered)

        # STEP 11: eta
        if config.calculate_eta:
            logger.info("STEP 11: Calculate eta values")
            _, kmax_df = get_eta(kapp_dfs_filtered, kmax_df)

        logger.success(f"p_total={p_label}: {len(kmax_df)} kmax entries computed")

        # Save individual kmax
        kmax_df.to_csv(output_dir / f"kmax_p{p_label}.csv", index=False)

        # Collect for merge: keep id cols + the kcat_app_max column renamed
        kmax_subset = kmax_df[id_cols + ["kcat_app_max"]].copy()
        kmax_subset = kmax_subset.rename(
            columns={"kcat_app_max": f"kcat_app_max_{p_label}"}
        )

        if merged is None:
            merged = kmax_subset
        else:
            merged = merged.merge(kmax_subset, on=id_cols, how="outer")

    # ==== Save merged sensitivity results ====
    logger.info("=" * 60)
    logger.info("Saving merged sensitivity results")
    out_path = output_dir / "kmax_sensitivity.csv"
    merged.to_csv(out_path, index=False)
    logger.success(f"Sensitivity results saved to: {out_path}")
    logger.info("=" * 60)

    return merged


def main():
    """Main entry point for the sensitivity analysis CLI."""
    parser = argparse.ArgumentParser(
        prog="run-sensitivity-analysis",
        description="Run p_total sensitivity analysis on the kapp pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python scripts/run_sensitivity_analysis.py configs/run_sensitivity_analysis/davidi_consensus.yaml

Output:
  Results are saved under results/run_sensitivity_analysis/{organism}_{folder_id}/results/
        """,
    )
    parser.add_argument(
        "config", type=Path, help="Path to sensitivity analysis YAML configuration file"
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true", help="Enable verbose (DEBUG) console output"
    )
    args = parser.parse_args()

    # Load sensitivity config
    try:
        sens_cfg = load_sensitivity_config(args.config)
    except FileNotFoundError as e:
        logger.error(f"Configuration file not found: {e}")
        sys.exit(1)
    except Exception as e:
        logger.error(f"Error loading sensitivity configuration: {e}")
        sys.exit(1)

    # Load base kapp pipeline config
    kapp_config_path = Path(sens_cfg["kapp_config"])
    try:
        logger.info(f"Loading base pipeline configuration from: {kapp_config_path}")
        config = load_config(kapp_config_path)
    except FileNotFoundError as e:
        logger.error(f"Base pipeline config not found: {e}")
        sys.exit(1)
    except Exception as e:
        logger.error(f"Error loading base pipeline config: {e}")
        sys.exit(1)

    # Build p_total sweep range
    p_min = sens_cfg.get("p_total_min", 0.26)
    p_max = sens_cfg.get("p_total_max", 0.31)
    p_step = sens_cfg.get("p_total_step", 0.01)
    p_total_values = list(np.round(np.arange(p_min, p_max + p_step / 2, p_step), 2))

    folder_id = sens_cfg.get("folder_id", config.folder_id)
    run_name = f"{config.organism}_{folder_id}"

    # Create output directories
    results_base = PROJ_ROOT / "results" / "run_sensitivity_analysis"
    run_root = results_base / run_name
    output_dir = ensure_dir_exists(run_root / "results")
    data_dir = ensure_dir_exists(run_root / "data")

    # Setup logging
    log_file = output_dir / f"log_{run_name}.log"
    setup_logging(log_file, run_name, verbose=args.verbose)

    # Log configuration summary
    logger.info("=" * 60)
    logger.info("SENSITIVITY ANALYSIS CONFIGURATION")
    logger.info("=" * 60)
    logger.info(f"Run date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    logger.info(f"Run name: {run_name}")
    logger.info(f"Organism: {config.organism}")
    logger.info(f"Base config: {kapp_config_path}")
    logger.info(f"Model: {config.model_path.name}")
    logger.info(f"Flux method: {config.flux_method}")
    logger.info(f"Solver: {config.solver}")
    logger.info(f"p_total sweep: {p_total_values}")
    logger.info(f"Output directory: {output_dir.relative_to(PROJ_ROOT)}")
    logger.info("=" * 60)

    # Run sensitivity analysis
    try:
        run_sensitivity_analysis(
            config=config,
            p_total_values=p_total_values,
            output_dir=output_dir,
            data_dir=data_dir,
            run_name=run_name,
        )
        logger.success("Sensitivity analysis completed successfully!")
        return 0
    except Exception as e:
        logger.exception(f"Sensitivity analysis failed: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
