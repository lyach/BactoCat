import pandas as pd
from pathlib import Path
from BactoCat import get_substrate_df, get_sequence_df

def run_kapp_pipeline(model_path: str, 
                      flux_method: str, 
                      carbon_uptake: list, 
                      oxygen_uptake: list, 
                      p_total: list,
                      substrate_df: str
                      sequence_df: str,
                      paxdb_dir: str):
    """
    Run the kapp pipeline.
    
    Parameters:
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
    
    Returns:
        kapp_dfs: dict
            Dictionary of kapp dataframes for each p_total value.
    """
    
    data_dir = Path(__file__).parent.parent / "data"
    results_dir = data_dir / "results"
    
    # Load the model
    try:
        model = cobra.io.read_sbml_model(model_path)
        model = model.copy() # Copy the model to avoid modifying the original
        print(f"Model loaded successfully from {model_path}")
    except Exception as e:
        raise ValueError(f"Model not found at {model_path}. Error: {e}")
    
    # ==== 1. Create GEM enzymes dataframe ====
    df_enzymes = create_gpr_dataframe(model)
    
    # Prints
    stats = analyze_model_gprs(model)
    print(f"\nModel Stats:")
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
    # TO MERGE LATER
    
    # ==== 4. Get sequence information ====
    if sequence_df:
        try:
            sequence_df = pd.read_csv(sequence_df)
        except:
            raise ValueError(f"Sequence dataframe not found at {sequence_df}"
    else:
        print(f"No sequence dataframe provided, retrieving sequences from UniProt.")
        sequence_df = get_sequence_df(model)
        
    # ==== 5. Get substrate information ====
    if substrate_df != None:
        try:
            substrate_df = pd.read_csv(substrate_df)
        except:
            raise ValueError(f"Substrate dataframe not found at {substrate_df}")
    else:
        substrate_df = get_substrate_df(model)
    
    # ==== 6. Create enzyme information dataframe ====
    enzymes_info_dfs = create_enzyme_info_dataframe(df_enzymes, fluxomics_df, substrates_df, sequence_df)
    print(f"Flux dataframe for condition 1 (carbon={carbon_uptake[0]}, oxygen={oxygen_uptake[0]}): {fluxomics_df['flux_cond1'].values[0]}")
    
    # === 7. Map proteomics information ====
    enzyme_protein_info_dfs = process_enzyme_protein_mapping(enzymes_info_dfs, paxdb_dir, p_total=p_total)
    
    # ==== 8. Calculate kapp for homomeric enzymes ====
    kapp_dfs = calculate_kapp_homomeric(enzyme_protein_info_dfs)
    
    # ==== 9. Filter values above physical threshold ====
    kapp_dfs_filtered = evaluate_kapp_homomeric(kapp_dfs)
    
    # ==== 10. Get kmax for homomeric enzymes ====
    kmax_dfs = get_kmax_homomeric(kapp_dfs_filtered)
    
    # ==== 11. Get eta values ====
    kapp_dfs_eta, kmax_dfs_eta_var = get_eta(kapp_dfs_filtered, kmax_dfs)
    
    # ==== 12. Save the results ====
    kmax_dfs_eta_var.to_csv(os.path.join(data_dir, "final", "kmax", "iml1515_homomeric_kmax_pFBA_variability.csv"), index=False)
    
    return kapp_dfs_eta, kmax_dfs_eta_var