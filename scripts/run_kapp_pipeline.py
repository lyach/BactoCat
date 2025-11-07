import pandas as pd
import cobra
import os
import sys
import yaml
from pathlib import Path

# Class to save console outputs
class Tee(object):
    def __init__(self, *files):
        self.files = files
    def write(self, obj):
        for f in self.files:
            f.write(str(obj))
            f.flush()
    def flush(self):
        for f in self.files:
            f.flush()

# Add project root to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

# Module imports
from src.enzyme_classifier import create_gpr_dataframe, analyze_model_gprs
# from src.gene_sequence_mapper import map_gem_genes_to_uniprot
# from src.substrate_mapper import get_substrate_df
from src.kapp_builder import (
    create_fluxomics_dataframe,
    create_enzyme_info_dataframe,
    process_enzyme_protein_mapping,
    calculate_kapp_homomeric,
    evaluate_kapp_homomeric,
    get_kmax_homomeric,
    get_eta
)

def run_kapp_pipeline(organism: str,
                      model_path: str,
                      flux_method: str, 
                      carbon_uptake: list, 
                      oxygen_uptake: list, 
                      p_total: list,
                      substrate_df: str = None,
                      sequence_df: str = None,
                      paxdb_path: str = None,
                      solver: str = "cplex",
                      output_dir: str = None,
                      data_dir: str = None
                      ):
    """
    Run the kapp pipeline.
    
    Parameters:
        organism: str
            Organism name.
        model_path: str
            Path to the model XML file.
        flux_method: str
            Method for the flux simulations: 'FBA' or 'pFBA'.
        carbon_uptake: list
            List of carbon uptake rates to test.
        oxygen_uptake: list
            List of oxygen uptake rates to test.
        p_total: list
            List of p_total values to test.
        substrate_df: str
            Path to substrate dataframe CSV file.
        sequence_df: str
            Path to sequence dataframe CSV file.
        paxdb_path: str
            Path to PaxDB proteomics data file.
        solver: str, optional
            Solver for the flux simulations. Default is 'cplex'.
        output_dir: str, optional
            Directory to save output files. Default is scripts/results/
    
    Returns:
        tuple: (kapp_dfs_eta, kmax_dfs_eta_var)
            - kapp_dfs_eta: dict of kapp dataframes with eta values
            - kmax_dfs_eta_var: DataFrame with kmax values and variance metrics
    """
    # Set solver
    try:
        cobra.Configuration().solver = solver
        print(f"Solver set to {solver}")
    except Exception as e:
        raise ValueError(f"Error setting solver: {e}")
    
    
    # Load the model
    try:
        model = cobra.io.read_sbml_model(model_path)
        model = model.copy()
        print(f"Model loaded successfully from {model_path}")
    except Exception as e:
        raise ValueError(f"Model not found at {model_path}. Error: {e}")
    
    # ==== 1. Create GEM enzymes dataframe ====
    df_enzymes = create_gpr_dataframe(model)
    
    # Prints
    stats = analyze_model_gprs(model)
    print("\nModel Stats:")
    print(f"Total reactions: {stats['total_reactions']}")
    print(f"Reactions with GPR: {stats['reactions_with_gpr']}")
    print(f"Total genes: {stats['total_genes']}")
    print(f"GPR cases: {stats['gpr_complexity']}")
    df_enzymes.head()
    
    # ==== 2. Get fluxomics simulations ====
    fluxomics_df = create_fluxomics_dataframe(flux_method=flux_method, GEM=model, 
                                         carbon_uptake=carbon_uptake, 
                                         oxygen_uptake=oxygen_uptake)
    
    # ==== 3. Get flux variability analysis ====
    from src.kapp_builder import create_FVA_dataframe, FVA_integration

    print("\nRunning Flux Variability Analysis (FVA)...")
    
    try:
        fva_df = create_FVA_dataframe(
            GEM_path=model_path,
            carbon_uptake=carbon_uptake,
            oxygen_uptake=oxygen_uptake,
            mu_fraction=0.9,
            solver=solver
        )
        print("FVA dataframe created successfully.")
    except Exception as e:
        raise RuntimeError(f"Error creating FVA dataframe: {e}") 
    
    print("\nIntegrating FVA results with fluxomics data...")
    try:
        filtered_fluxomics_df, violations_df = FVA_integration(fluxomics_df, fva_df, filter=True)
        fluxomics_df = filtered_fluxomics_df.copy()

        # Save outputs for reference
        fva_df.to_csv(output_dir / "FVA_bounds.csv", index=False)
        filtered_fluxomics_df.to_csv(output_dir / "fluxomics_filtered.csv", index=False)
        violations_df.to_csv(output_dir / "FVA_violations.csv", index=False)

        print(f"FVA integration complete. Filtered fluxomics: {filtered_fluxomics_df.shape[0]} rows")
        print(f"Violations detected: {violations_df.shape[0]} rows")
    except Exception as e:
        raise RuntimeError(f"Error during FVA integration: {e}")

    
    # ==== 4. Get sequence information ====
    from src.gene_sequence_mapper import map_organism_to_uniprot
    print("\n==== 4. Loading sequence information ====")
    if sequence_df:
        try:
            sequence_df_loaded = pd.read_csv(sequence_df)
            print(f"Sequence dataframe loaded: {len(sequence_df_loaded)} rows")
        except Exception as e:
            raise ValueError(f"Sequence dataframe not found at {sequence_df}. Error: {e}")
    else: 
        print(f"No sequence dataframe provided, retrieving sequences from UniProt.")
        sequence_df_loaded = map_organism_to_uniprot(organism)

    sequence_df_loaded.to_csv(data_dir / "sequence_df.csv", index=False)

        
    # ==== 5. Get substrate information ====
    print("\n==== 5. Loading substrate information ====")
    try:
        substrate_df_loaded = pd.read_csv(substrate_df)
        print(f"Substrate dataframe loaded: {len(substrate_df_loaded)} rows")
    except Exception as e:
        raise ValueError(f"Substrate dataframe not found at {substrate_df}. Error: {e}")
    
    # if substrate_df:
    #     try:
    #         substrate_df = pd.read_csv(substrate_df)
    #     except:
    #         raise ValueError(f"Substrate dataframe not found at {substrate_df}")
    # else:
    #     substrate_df = get_substrate_df(model)

    
    # ==== 6. Create enzyme information dataframe ====
    print("\n==== 6. Creating enzyme information dataframe ====")
    enzymes_info_dfs = create_enzyme_info_dataframe(df_enzymes, fluxomics_df, substrate_df_loaded, sequence_df_loaded)
    
    # === 7. Map proteomics information ====
    print("\n==== 7. Mapping proteomics information ====")
    enzyme_protein_info_dfs = process_enzyme_protein_mapping(enzymes_info_dfs, paxdb_path, p_total=p_total)
    
    # ==== 8. Calculate kapp for homomeric enzymes ====
    print("\n==== 8. Calculating kapp for homomeric enzymes ====")
    kapp_dfs = calculate_kapp_homomeric(enzyme_protein_info_dfs)
    
    # ==== 9. Filter values above physical threshold ====
    kapp_dfs_filtered = evaluate_kapp_homomeric(kapp_dfs)
    
    # ==== 10. Get kmax for homomeric enzymes ====
    print("\n==== 10. Getting kmax for homomeric enzymes ====")
    kmax_dfs = get_kmax_homomeric(kapp_dfs_filtered)
    
    # ==== 11. Get eta values ====
    print("\n==== 11. Calculating eta values ====")
    kapp_dfs_eta, kmax_dfs_eta_var = get_eta(kapp_dfs_filtered, kmax_dfs)
    
    # ==== 12. Save the results ====
    print("\n==== 12. Saving results ====")
    output_file = output_dir / f"iml1515_homomeric_kmax_{flux_method}_variability.csv"
    kmax_dfs_eta_var.to_csv(output_file, index=False)
    print(f"Results saved to: {output_file}")
    
    print(f"\n{'='*60}")
    print("Pipeline completed!")
    print(f"{'='*60}")
    
    return kapp_dfs_eta, kmax_dfs_eta_var


def main():
    """
    Main function to run the pipeline from command line with config file.
    """
    import argparse
    
    # Set up argument parser
    parser = argparse.ArgumentParser(
        description='Run kapp pipeline for a metabolic model',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
        Example usage:
        python run_kapp_pipeline.py config.yaml
        python run_kapp_pipeline.py config.yaml --output-dir ./custom_results
        """
    )
    
    parser.add_argument(
        'config',
        type=str,
        help='Path to YAML configuration file'
    )
    
    parser.add_argument(
        '--output-dir',
        type=str,
        default=None,
        help='Override output directory from config file'
    )
    
    args = parser.parse_args()
    
    # Load configuration file
    print(f"Loading configuration from: {args.config}")
    try:
        with open(args.config, 'r', encoding='utf-8') as f:
            config = yaml.safe_load(f)
    except FileNotFoundError:
        print(f"Error: Configuration file not found: {args.config}")
        sys.exit(1)
    except yaml.YAMLError as e:
        print(f"Error parsing YAML configuration: {e}")
        sys.exit(1)
        
    script_dir = Path(__file__).parent
    
    # Extract parameters from config
    try:
        organism = config['organism']
        model_path = script_dir / config['model_path']
        flux_method = config['flux_method']
        carbon_uptake = config['carbon_uptake']
        oxygen_uptake = config['oxygen_uptake']
        p_total = config['p_total']
        paxdb_path = script_dir / config['paxdb_path']
        solver = config.get('solver', 'cplex')  # Default to 'cplex' if not specified
        output_dir_str = args.output_dir or config.get('output_dir', None)
        # Optional parameters
        raw_substrate_df = config.get('substrate_df')
        raw_sequence_df = config.get('sequence_df')
        substrate_df = script_dir / raw_substrate_df if raw_substrate_df else None
        sequence_df = script_dir / raw_sequence_df if raw_sequence_df else None
    except KeyError as e:
        print(f"Error: Missing required parameter in config file: {e}")
        sys.exit(1)
        
    # Resolve output directory
    if output_dir_str is None:
        # Default root directory is under scripts/results/
        run_root = script_dir / f"{organism}"
    else:
        # Use user-provided path as the root
        run_root_path = Path(output_dir_str)
        if not run_root_path.is_absolute():
            run_root = script_dir / run_root_path
        else:
            run_root = run_root_path
    
    # Define and create the structured subdirectories
    output_dir = run_root / "results" # Where config, logs, and final outputs go
    data_dir = run_root / "data"     # Where intermediate data (sequence_df, FVA) goes

    # Create directories if they don't exist
    output_dir.mkdir(parents=True, exist_ok=True)
    data_dir.mkdir(parents=True, exist_ok=True)
    print(f"Output directory set to: {output_dir.relative_to(script_dir).as_posix()}")
    print(f"Data directory set to: {data_dir.relative_to(script_dir).as_posix()}")

    # Helper lambda to safely format the path for logging
    path_to_log = lambda p: p.relative_to(script_dir).as_posix() if p else 'Auto-generated'
    
    # Display and save configuration
    config_lines = [
        f"\n{'='*60}",
        "KAPP PIPELINE CONFIGURATION",
        f"{'='*60}",
        f" Organism: {organism}",
        f" Model: {model_path.relative_to(script_dir).as_posix()}",
        f" Flux method: {flux_method}",
        f" Solver: {solver}",
        f" Carbon uptake rates: {carbon_uptake}",
        f" Oxygen uptake rates: {oxygen_uptake}",
        f" P_total values: {p_total}",
        f" Substrate data: {path_to_log(substrate_df)}",
        f" Sequence data: {path_to_log(sequence_df)}",
        f" PaxDB data: {paxdb_path.relative_to(script_dir).as_posix()}",
        f" Output directory: {output_dir.relative_to(script_dir).as_posix()}",
        f"{'='*60}\n"
    ]
    config_text = "\n".join(config_lines)
    print(config_text)
    
    config_filepath = output_dir / f"kmax_homomeric_{organism}.txt"
    with open(config_filepath, 'w', encoding='utf-8') as f:
        f.write(config_text)
    
    # Redirect prints to the log file and console
    original_stdout = sys.stdout
    with open(config_filepath, 'a', encoding='utf-8') as log_file:
        sys.stdout = Tee(sys.stdout, log_file)
        
        try:
            print(f"\nConfiguration saved to {config_filepath.relative_to(script_dir).as_posix()}")

            # Run the pipeline
            kapp_results, kmax_results = run_kapp_pipeline(
                organism=organism,
                model_path=model_path,
                flux_method=flux_method,
                carbon_uptake=carbon_uptake,
                oxygen_uptake=oxygen_uptake,
                p_total=p_total,
                substrate_df=substrate_df,
                sequence_df=sequence_df,
                paxdb_path=paxdb_path,
                solver=solver,
                output_dir=output_dir,
                data_dir=data_dir
            )
            
            print("\nPipeline execution completed successfully!")
                    
        except Exception as e:
            print(f"\nError running pipeline: {e}")
            import traceback
            traceback.print_exc()
            sys.exit(1)
        finally:
            # Restore stdout
            sys.stdout = original_stdout


if __name__ == "__main__":
    main()