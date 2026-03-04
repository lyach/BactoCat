"""
BactoCat Kapp Pipeline Runner.

This script runs the kapp pipeline for building enzyme kinetic datasets
from genome-scale metabolic models.

Usage:
    run-kapp-pipeline config.yaml
    python -m scripts.run_kapp_pipeline config.yaml
"""

# Command line interface
import argparse
import sys
from datetime import datetime
from pathlib import Path
import random

# Scientific libraries
import cobra
import pandas as pd
from loguru import logger

# Module imports
from src.config import PipelineConfig, load_config, PROJ_ROOT, ensure_dir_exists
from src.enzyme_classifier import create_gpr_dataframe, analyze_model_gprs
from src.gene_sequence_mapper import map_organism_to_uniprot
from src.substrate_mapper import get_substrate_df
from src.kapp_builder import (
    create_fluxomics_dataframe,
    create_enzyme_info_dataframe,
    process_enzyme_protein_mapping,
    create_FVA_dataframe,
    FVA_integration,
    calculate_kapp_homomeric,
    evaluate_kapp_homomeric,
    get_kmax_homomeric,
    get_eta,
)


def setup_logging(log_file: Path, run_name: str) -> None:
    """
    Configure loguru for console and file output.
    
    Parameters
    ----------
    log_file : Path
        Path to the log file
    run_name : str
        Name of the run (for log formatting)
    """
    # Remove default handler
    logger.remove()
    
    # Console handler with colored output
    logger.add(
        sys.stderr,
        format="<green>{time:HH:mm:ss}</green> | <level>{level: <8}</level> | <level>{message}</level>",
        level="INFO",
        colorize=True,
    )
    
    # File handler with full details
    logger.add(
        log_file,
        format="{time:YYYY-MM-DD HH:mm:ss} | {level: <8} | {message}",
        level="DEBUG",
        rotation="10 MB",
        retention="30 days",
    )
    
    logger.info(f"Logging initialized for run: {run_name}")


def run_kapp_pipeline(
    config: PipelineConfig,
    output_dir: Path,
    data_dir: Path,
    run_name: str,
) -> tuple[dict, pd.DataFrame]:
    """
    Run the kapp pipeline.
    
    Parameters
    ----------
    config : PipelineConfig
        Validated pipeline configuration
    output_dir : Path
        Directory for output files
    data_dir : Path
        Directory for intermediate data files
    run_name : str
        Name for this run (used in file names)
        
    Returns
    -------
    tuple
        (kapp_dfs_eta, kmax_dfs_eta_var)
        - kapp_dfs_eta: dict of kapp DataFrames with eta values
        - kmax_dfs_eta_var: DataFrame with kmax values and variance metrics
    """
    # Set solver
    try:
        cobra.Configuration().solver = config.solver
        logger.info(f"Solver set to {config.solver}")
    except Exception as e:
        raise ValueError(f"Error setting solver: {e}")
    
    # Load the model
    try:
        model = cobra.io.read_sbml_model(str(config.model_path))
        model = model.copy()
        logger.info(f"Model loaded successfully from {config.model_path.name}")
        
        # Log solver tolerance settings
        try:
            solver_instance = model.solver
            if hasattr(solver_instance, 'configuration') and hasattr(solver_instance.configuration, 'tolerances'):
                tolerances = solver_instance.configuration.tolerances
                logger.debug(f"Feasibility Tolerance: {tolerances.feasibility}")
                logger.debug(f"Optimality Tolerance: {tolerances.optimality}")
                logger.debug(f"Integrality Tolerance: {tolerances.integrality}")
        except Exception:
            logger.debug("Solver tolerance information not available")
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
    logger.debug(f"GPR cases: {stats['gpr_complexity']}")
    
    # ==== STEP 2: Run fluxomics simulations ====
    logger.info("=" * 50)
    logger.info(f"STEP 2: Run {config.flux_method} fluxomics simulations")
    
    # Load medium_df if provided
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
        carbon_uptake=config.carbon_uptake,
        oxygen_uptake=config.oxygen_uptake,
        carbon_exchange_rxn=config.carbon_exchange_rxn,
        oxygen_exchange_rxn=config.oxygen_exchange_rxn,
        medium_df=medium_df_loaded,
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
                carbon_uptake=config.carbon_uptake,
                oxygen_uptake=config.oxygen_uptake,
                mu_fraction=config.mu_fraction,
                solver=config.solver,
                carbon_exchange_rxn=config.carbon_exchange_rxn,
                oxygen_exchange_rxn=config.oxygen_exchange_rxn,
                medium_df=medium_df_loaded,
            )
            logger.info("FVA dataframe created successfully")

            logger.info("Integrating FVA results with fluxomics data")
            filtered_fluxomics_df, violations_df = FVA_integration(fluxomics_df, fva_df, filter=True)
            fva_df.to_csv(output_dir / f"FVA_bounds_{run_name}.csv", index=False)
            filtered_fluxomics_df.to_csv(output_dir / f"fluxomics_filtered_{run_name}.csv", index=False)
            violations_df.to_csv(output_dir / f"FVA_violations_{run_name}.csv", index=False)

            logger.info(f"FVA integration complete. Filtered fluxomics: {filtered_fluxomics_df.shape[0]} rows")
            logger.info(f"Violations detected: {violations_df.shape[0]} rows")
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
        logger.info(f"Sequence dataframe created: {len(sequence_df_loaded)} rows at {seq_output.name}")
    
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
        logger.info(f"Substrate dataframe created: {len(substrate_df_loaded)} rows at {sub_output.name}")
    
    # ==== STEP 6: Create enzyme information dataframe ====
    logger.info("=" * 50)
    logger.info("STEP 6: Create enzyme information dataframe")
    enzymes_info_dfs = create_enzyme_info_dataframe(
        df_enzymes, fluxomics_df, substrate_df_loaded, sequence_df_loaded, run_fva=config.run_fva
    )
    
    # ==== STEP 7: Map proteomics information ====
    logger.info("=" * 50)
    logger.info("STEP 7: Map proteomics information")
    enzyme_protein_info_dfs = process_enzyme_protein_mapping(
        enzymes_info_dfs, str(config.paxdb_path), p_total=config.p_total
    )
    
    # ==== STEP 8: Calculate kapp for homomeric enzymes ====
    logger.info("=" * 50)
    logger.info("STEP 8: Calculate kapp for homomeric enzymes")
    kapp_dfs = calculate_kapp_homomeric(enzyme_protein_info_dfs)
    
    # ==== STEP 9: Filter values above physical threshold ====
    logger.info("=" * 50)
    logger.info("STEP 9: Filter values above physical threshold")
    kapp_dfs_filtered = evaluate_kapp_homomeric(kapp_dfs)
    
    # ==== STEP 10: Get kmax for homomeric enzymes ====
    logger.info("=" * 50)
    logger.info("STEP 10: Get kmax for homomeric enzymes")
    kmax_df = get_kmax_homomeric(kapp_dfs_filtered)
    
    # ==== STEP 11: Calculate eta values ====
    if config.calculate_eta:
        logger.info("=" * 50)
        logger.info("STEP 11: Calculate eta values")
        kapp_dfs_eta, kmax_df_out = get_eta(kapp_dfs_filtered, kmax_df)
    else:
        kapp_dfs_eta = kapp_dfs_filtered
        kmax_df_out = kmax_df
    
    # ==== STEP 12: Save results ====
    logger.info("=" * 50)
    logger.info("STEP 12: Save results")
    output_file = output_dir / f"kmax.csv"
    kmax_df_out.to_csv(output_file, index=False)
    logger.success(f"Results saved to: {output_file}")
    
    logger.info("=" * 60)
    logger.success("Pipeline completed!")
    logger.info("=" * 60)
    
    return kmax_df_out


def main():
    """Main entry point for the kapp pipeline CLI."""
    parser = argparse.ArgumentParser(
        prog="run-kapp-pipeline",
        description="Run kapp pipeline for building enzyme kinetic datasets from metabolic models",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  run-kapp-pipeline configs/run_kapp_pipeline/ecoli_homomeric.yaml
  python -m scripts.run_kapp_pipeline configs/run_kapp_pipeline/ecoli_homomeric.yaml

Output:
  Results are saved to results/run_kapp_pipeline/{organism}_{date}_{id}/
  with subdirectories /data (intermediate files) and /results (final outputs)
        """,
    )
    
    parser.add_argument(
        "config",
        type=Path,
        help="Path to YAML configuration file",
    )
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Enable verbose (DEBUG) console output",
    )
    
    args = parser.parse_args()
    
    # Load and validate configuration
    try:
        logger.info(f"Loading configuration from: {args.config}")
        config = load_config(args.config)
    except FileNotFoundError as e:
        logger.error(f"Configuration file not found: {e}")
        sys.exit(1)
    except Exception as e:
        logger.error(f"Error loading configuration: {e}")
        sys.exit(1)
    
    # Set folder name
    run_name = f"{config.organism}_{config.folder_id}"
    
    # Create output directories
    results_base = PROJ_ROOT / "results" / "run_kapp_pipeline"
    run_root = results_base / run_name
    output_dir = ensure_dir_exists(run_root / "results")
    data_dir = ensure_dir_exists(run_root / "data")
    
    # Setup logging with file output
    log_file = output_dir / f"log_{run_name}.log"
    setup_logging(log_file, run_name)
    
    # Adjust console log level if verbose
    if args.verbose:
        logger.remove()
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
    
    # Log configuration summary
    logger.info("=" * 60)
    logger.info("KAPP PIPELINE CONFIGURATION")
    logger.info("=" * 60)
    logger.info(f"Run name: {run_name}")
    logger.info(f"Organism: {config.organism}")
    logger.info(f"Model: {config.model_path.name}")
    logger.info(f"Flux method: {config.flux_method}")
    logger.info(f"Solver: {config.solver}")
    if config.medium_df is not None:
        logger.info(f"Medium dataframe: {config.medium_df.name}")
    else:
        logger.info(f"Carbon uptake rates: {config.carbon_uptake}")
        logger.info(f"Oxygen uptake rates: {config.oxygen_uptake}")
        logger.info(f"Carbon exchange rxn: {config.carbon_exchange_rxn}")
        logger.info(f"Oxygen exchange rxn: {config.oxygen_exchange_rxn}")
    logger.info(f"P_total: {config.p_total}")
    logger.info(f"Substrate data: {config.substrate_df.name if config.substrate_df else 'Auto-generated'}")
    logger.info(f"Sequence data: {config.sequence_df.name if config.sequence_df else 'Auto-generated'}")
    logger.info(f"PaxDB data: {config.paxdb_path.name}")
    logger.info(f"Output directory: {output_dir.relative_to(PROJ_ROOT)}")
    logger.info(f"Data directory: {data_dir.relative_to(PROJ_ROOT)}")
    logger.info("=" * 60)
    
    # Run the pipeline
    try:
        kmax_result = run_kapp_pipeline(
            config=config,
            output_dir=output_dir,
            data_dir=data_dir,
            run_name=run_name,
        )
        logger.success("Pipeline execution completed successfully!")
        return 0
        
    except Exception as e:
        logger.exception(f"Pipeline failed with error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
