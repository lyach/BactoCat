"""
Module Description: kapp_builder.py

Purpose: 
Script for building the kapp dataframe for homomeric and heteromeric enzymes.

Overview: 
This module provides functions to:
- Calculate apparent kcat (kapp) for homomeric enzymes
- Calculate specific activity (SA_app) for heteromeric enzyme complexes
"""

import cobra
import pandas as pd
import numpy as np
from cobra import flux_analysis
from itertools import product
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from cobra.io import read_sbml_model
from src.FVA_analysis.utils import cobra_to_fva_problem
import json

def create_fluxomics_dataframe(flux_method: str, GEM: cobra.Model, 
                              carbon_uptake: list, oxygen_uptake: list):
    """
    Create a dataframe with FBA or pFBA fluxomics for all combinations of carbon and oxygen uptake.
    
    Parameters:
        flux_method: str
            Method for the flux simulations: 'FBA' or 'pFBA'
        GEM: cobra.Model
            The GEM model to perform flux analysis on
        carbon_uptake: list
            List of carbon uptake rates to test
        oxygen_uptake: list
            List of oxygen uptake rates to test
    
    Returns:
        pd.DataFrame: DataFrame with columns: 'rxn_id', 'flux_cond1', 'flux_cond2', ..., 'flux_condN'
    """
    
    # Create all combinations of carbon and oxygen uptake rates
    uptake_combinations = list(product(carbon_uptake, oxygen_uptake))
    
    # Initialize results dictionary
    flux_results = {}
    
    # Get all reaction IDs for the dataframe
    rxn_ids = [rxn.id for rxn in GEM.reactions]
    
    # Process each combination
    for i, (carbon_rate, oxygen_rate) in enumerate(uptake_combinations, 1):
        print(f"Processing condition {i}: Carbon={carbon_rate}, Oxygen={oxygen_rate}")
        
        # Create a copy of the model to avoid modifying the original
        model_copy = GEM.copy()
        
        # Set carbon uptake rate
        try:
            carbon_rxn = model_copy.reactions.get_by_id('EX_glc__D_e')  # TO DO - add input param to specify the carbon uptake reaction
            carbon_rxn.lower_bound = -abs(carbon_rate)  # Negative for uptake
        except KeyError:
            print("Warning: Carbon uptake reaction not found. Skipping carbon constraint.")
        
        # Set oxygen uptake rate
        try:
            oxygen_rxn = model_copy.reactions.get_by_id('EX_o2_e')  # TO DO - add input param to specify the oxygen uptake reaction
            oxygen_rxn.lower_bound = -abs(oxygen_rate)  # Negative for uptake
        except KeyError:
            print("Warning: Oxygen uptake reaction not found. Skipping oxygen constraint.")
        
        # Run FBA or pFBA
        if flux_method == 'FBA':
            solution = model_copy.optimize()
        elif flux_method == 'pFBA':
            solution = flux_analysis.pfba(model_copy)
        else:
            raise ValueError(f"Invalid method '{flux_method}'. Must be 'FBA' or 'pFBA'.")
        
        # Store results
        if solution.status == 'optimal':
            flux_results[f'flux_cond{i}'] = [solution.fluxes[rxn_id] for rxn_id in rxn_ids]
            print(f"Condition {i} completed successfully")
        else:
            print(f"Warning: Condition {i} optimization failed with status: {solution.status}")
            flux_results[f'flux_cond{i}'] = [0.0] * len(rxn_ids)  # Fill with zeros for failed optimization
    
    # Create the output dataframe
    fluxomics_df = pd.DataFrame({'rxn_id': rxn_ids})
    
    # Add flux columns for each condition
    for condition, fluxes in flux_results.items():
        fluxomics_df[condition] = fluxes
    
    print(f"Fluxomics dataframe created with {len(uptake_combinations)} conditions")
    return fluxomics_df

def create_FVA_dataframe(GEM_path: str, 
                                        carbon_uptake: list, 
                                        oxygen_uptake: list,
                                        mu_fraction: float = 0.9,
                                        solver: str = 'cplex'):
    """
    Run FVA for all combinations of carbon and oxygen uptake rates, 
    matching the structure of create_fluxomics_dataframe().
    
    Parameters
    ----------
    GEM_path : str
        Path to the SBML model file (XML).
    carbon_uptake : list
        List of carbon uptake rates to test.
    oxygen_uptake : list
        List of oxygen uptake rates to test.
    mu_fraction : float, optional
        Fraction of optimal growth rate for FVA (default = 0.9).
    
    Returns
    -------
    pd.DataFrame
        Combined FVA dataframe with columns:
        ['rxn_id', 'FVA_lower_cond1', 'FVA_upper_cond1', ..., 'FVA_lower_condN', 'FVA_upper_condN']
    """
    
    # Conditionally import the correct FVA solver
    if solver.lower() == 'cplex':
        from src.FVA_analysis.fvfa_cplex import fva_solve_faster
    elif solver.lower() == 'gurobi':
        from src.FVA_analysis.fvfa import fva_solve_faster
    else:
        raise ValueError(f"Solver '{solver}' is not supported for FVA. Please use 'cplex' or 'gurobi'.")

    # Create all combinations of carbon and oxygen uptake rates
    uptake_combinations = list(product(carbon_uptake, oxygen_uptake))
    
    # Load base model once
    base_model = read_sbml_model(GEM_path)
    rxn_ids = [rxn.id for rxn in base_model.reactions]
    
    # Initialize dictionaries for lower/upper bounds
    FVA_lower_results = {}
    FVA_upper_results = {}
    
    for i, (carbon_rate, oxygen_rate) in enumerate(uptake_combinations, 1):
        print(f"Running FVA condition {i}: Carbon={carbon_rate}, Oxygen={oxygen_rate}")
        
        # Copy model to avoid media interference
        model_copy = base_model.copy()
        
        # Set medium
        try:
            model_copy.reactions.EX_glc__D_e.lower_bound = -abs(carbon_rate)
        except KeyError:
            print("Warning: Carbon uptake reaction not found. Skipping.")
        try:
            model_copy.reactions.EX_o2_e.lower_bound = -abs(oxygen_rate)
        except KeyError:
            print("Warning: Oxygen uptake reaction not found. Skipping.")
        
        # Optimize and get optimal mu
        solution = model_copy.optimize()
        if solution.status != 'optimal':
            print(f"Warning: optimization failed at condition {i} with status:({solution.status}), filling with NaNs.")
            FVA_lower_results[f'FVA_lower_cond{i}'] = [float('nan')] * len(rxn_ids)
            FVA_upper_results[f'FVA_upper_cond{i}'] = [float('nan')] * len(rxn_ids)
            continue
        
        # Build FVA problem
        problem = cobra_to_fva_problem(model_copy, mu=mu_fraction)
        
        # Run FVA
        fva_results = fva_solve_faster(problem)
        rxn_ids = [rxn.id for rxn in model_copy.reactions]
        fva_df = pd.DataFrame({
        'rxn_id': [rxn.id for rxn in model_copy.reactions],
        'FVA_lower': fva_results.lower_bound,
        'FVA_upper': fva_results.upper_bound
    }) 
        fva_df = fva_df.set_index('rxn_id').reindex(rxn_ids).reset_index()
        
        # Store FVA lower/upper bounds
        FVA_lower_results[f'FVA_lower_cond{i}'] = fva_df['FVA_lower'].values
        FVA_upper_results[f'FVA_upper_cond{i}'] = fva_df['FVA_upper'].values
        
        print(f"Condition {i} completed successfully.")
    
    # Build the output dataframe
    fva_combined = pd.DataFrame({'rxn_id': rxn_ids})
    
    for col_name, values in FVA_lower_results.items():
        fva_combined[col_name] = values
    for col_name, values in FVA_upper_results.items():
        fva_combined[col_name] = values
    
    print(f"FVA dataframe created with {len(uptake_combinations)} conditions.")
    return fva_combined


def FVA_integration(fluxomics_df: pd.DataFrame, fva_df: pd.DataFrame, filter: bool = False):
    """
    Check fluxes against FVA bounds for all conditions and optionally filter out reactions with violations.
    
    Returns:
        filtered_fluxomics_df, violations_df
    """
    merged_df = fluxomics_df.merge(fva_df, on='rxn_id', how='left')
    
    # Identify flux columns
    flux_cols = [col for col in merged_df.columns if col.startswith('flux_cond')]
    
    violations = []
    
    for col in flux_cols:
        # Map flux column to corresponding FVA lower/upper columns
        cond_num = col.split('flux_cond')[-1]
        lower_col = f'FVA_lower_cond{cond_num}'
        upper_col = f'FVA_upper_cond{cond_num}'
        
        below_mask = merged_df[col] < merged_df[lower_col]
        above_mask = merged_df[col] > merged_df[upper_col]
        
        for idx in merged_df[below_mask].index:
            violations.append({
                'rxn_id': merged_df.at[idx, 'rxn_id'],
                'condition': col,
                'flux': merged_df.at[idx, col],
                'FVA_lower': merged_df.at[idx, lower_col],
                'FVA_upper': merged_df.at[idx, upper_col],
                'violation_type': 'below_min'
            })
        for idx in merged_df[above_mask].index:
            violations.append({
                'rxn_id': merged_df.at[idx, 'rxn_id'],
                'condition': col,
                'flux': merged_df.at[idx, col],
                'FVA_lower': merged_df.at[idx, lower_col],
                'FVA_upper': merged_df.at[idx, upper_col],
                'violation_type': 'above_max'
            })
    
    violations_df = pd.DataFrame(violations)
    print(f"Detected {len(violations_df)} violations of FVA bounds.")
    
    if filter and not violations_df.empty:
        violating_rxns = violations_df['rxn_id'].unique()
        before = len(merged_df)
        merged_df = merged_df[~merged_df['rxn_id'].isin(violating_rxns)].copy()
        after = len(merged_df)
        print(f"Filtered out {before - after} reactions with violations.")
    
    filtered_fluxomics_df = merged_df.copy()
    return filtered_fluxomics_df, violations_df

def load_dataframe_if_path(data_input):
    """
    Auxiliary function to load dataframe from CSV path or return existing dataframe.
    
    Parameters:
        data_input: str or pd.DataFrame
            Either a file path to CSV or an existing DataFrame
    
    Returns:
        pd.DataFrame: The loaded or existing DataFrame
    """
    if isinstance(data_input, str):
        return pd.read_csv(data_input)
    elif isinstance(data_input, pd.DataFrame):
        return data_input
    else:
        raise ValueError("Input must be either a file path (str) or a pandas DataFrame")


def create_enzyme_info_dataframe(enzymes_df, fluxomics_df, substrates_df, sequence_df):
    """
    Create enzyme info dataframes for each flux condition by merging enzyme, flux, substrate, and sequence data.
    
    Parameters:
        enzymes_df: str or pd.DataFrame
            Enzymes dataframe or path to CSV file
        fluxomics_df: str or pd.DataFrame
            Fluxomics dataframe or path to CSV file
        substrates_df: str or pd.DataFrame
            Substrates dataframe or path to CSV file
        sequence_df: str or pd.DataFrame
            Sequence dataframe or path to CSV file
    
    Returns:
        enzyme_info_dfs: dict
            Dictionary with keys as condition names and values as processed DataFrames
    """
    
    # Load dataframes if paths are provided
    enzymes_df = load_dataframe_if_path(enzymes_df)
    fluxomics_df = load_dataframe_if_path(fluxomics_df)
    substrates_df = load_dataframe_if_path(substrates_df)
    sequence_df = load_dataframe_if_path(sequence_df)
    
    # Initialize output dictionary
    enzyme_info_dfs = {}
    
    # Get all flux condition columns
    flux_columns = [col for col in fluxomics_df.columns if col.startswith('flux_cond')]
    
    print(f"Processing {len(flux_columns)} flux conditions...")
    
    # Process each flux condition
    for flux_col in flux_columns:
        print(f"Processing {flux_col}...")
        print(f"Rows before filtering: {len(enzymes_df)}")
        
        # Create a copy of enzymes_df for this condition
        condition_df = enzymes_df.copy()
        cond_num = flux_col.replace("flux_cond", "")
        lower_col = f"FVA_lower_cond{cond_num}"
        upper_col = f"FVA_upper_cond{cond_num}"


        # == Fluxomics info ==
        # Merge with specific flux condition
        if lower_col not in fluxomics_df.columns or upper_col not in fluxomics_df.columns:
            print(f"Warning: Missing FVA columns for {flux_col} — skipping FVA merge.")
            flux_subset = fluxomics_df[['rxn_id', flux_col]].copy()
            flux_subset.rename(columns={flux_col: 'flux_value'}, inplace=True)
        else:
            flux_subset = fluxomics_df[['rxn_id', flux_col, lower_col, upper_col]].copy()
            flux_subset.rename(columns={
                flux_col: 'flux_value',
                lower_col: 'FVA_lower',
                upper_col: 'FVA_upper'
        }, inplace=True)
        
        
        condition_df = pd.merge(condition_df, flux_subset, left_on="rxn", right_on="rxn_id", how="left")
        
        # == Substrate info ==
        # Clean substrates df - kinGEMs specific
        substrates_clean = substrates_df[['Reaction', 'SMILES', 'Direction']].copy()
        
        # Merge substrate info
        condition_df = pd.merge(condition_df, substrates_clean, left_on="rxn", right_on="Reaction", how="left")
        
        # == Sequence info ==
        # Merge sequence info
        condition_df = pd.merge(condition_df, sequence_df, left_on="gene", right_on="model_gene_id", how="left")
        
        # == Data cleaning and filtering ==
        # Drop rows with wrong direction-flux
        condition_df = condition_df[
            ((condition_df['Direction'] == 'forward') & (condition_df['flux_value'] >= 0)) |
            ((condition_df['Direction'] == 'reverse') & (condition_df['flux_value'] <= 0))
        ]
        
        # Drop rows with balancing species / cofactors (unlikely to be substrates)
        cofactors = {
            # Common cofactors
            'O',           # water
            'O=O',         # molecular oxygen  
            '[H]',         # hydrogen
            'OO',          # hydrogen peroxide
            '[O]',         # atomic oxygen
            '[OH-]',       # hydroxide
            '[H+]',        # proton
            'N',           # nitrogen (sometimes used)
            'P',           # phosphorus
            'S',           # sulfur
            # Simple cofactors
            #'C(=O)O',      # formic acid
            #'CO',          # methanol
            #'CCO',         # ethanol
            #'CC(=O)O',     # acetic acid
            #'C',           # methane
            #'CC',          # ethane
            #'CCC',         # propane
            #'N',           # ammonia (as N)
            'NN',          # hydrazine
            'C=O',         # formaldehyde
            'CC=O',        # acetaldehyde
            'O=C=O',       # carbon dioxide
            '[NH3+]',      # ammonium
            '[Na+]',       # sodium
            '[Cl-]',       # chloride
            '[K+]',        # potassium
            '[Mg+2]',      # magnesium
            '[Ca+2]',      # calcium
        }
        condition_df = condition_df[~condition_df['SMILES'].isin(cofactors)]
        
        # Remove rows with 0 flux 
        condition_df = condition_df[condition_df['flux_value'] != 0]
        
        # Store the processed dataframe
        enzyme_info_dfs[flux_col] = condition_df
        
        print(f"Completed {flux_col}: {len(condition_df)} rows after filtering")
    
    print(f"Created enzyme info dataframes for {len(enzyme_info_dfs)} conditions")
    return enzyme_info_dfs


def calculate_molecular_weight(sequence: str) -> float:
    """
    Calculate molecular weight of a protein sequence using BioPython.
    
    Parameters
        sequence : str
            Amino acid sequence
    Returns
        float
            Molecular weight in g/mol (Daltons)
    """
    try:
        if pd.isna(sequence) or not sequence.strip():
            return float('nan')
        analysed_seq = ProteinAnalysis(sequence.strip())
        return analysed_seq.molecular_weight()
    except Exception:
        return float('nan')


def map_paxdb_to_gene(paxdb_df: pd.DataFrame, df_enzymes: pd.DataFrame, p_total: float) -> pd.DataFrame:
    """
    Map PaxDB abundances to enzymes by gene ID and calculate protein concentrations.

    Parameters
    paxdb_df : pd.DataFrame
        PaxDB dataframe with 'string_external_id' and 'abundance' columns
    df_enzymes : pd.DataFrame
        Enzyme dataframe with 'gene' and 'sequence' columns
    p_total : float
        Total protein content in g/gDCW
    
    Returns
    pd.DataFrame
        Copy of df_enzymes with new columns:
        - 'protein_ppm' (float, NaN if no match)
        - 'protein_mmol_gdcw' (float, protein concentration in mmol/gDCW)
    """
    # Work on copies
    pax = paxdb_df.copy()
    enz = df_enzymes.copy()

    # Ensure required columns exist
    required_pax_cols = {"string_external_id", "abundance"}
    if not required_pax_cols.issubset(pax.columns):
        missing = required_pax_cols - set(pax.columns)
        raise KeyError(f"paxdb_df missing columns: {missing}")
    
    required_enz_cols = {"gene", "sequence"}
    if not required_enz_cols.issubset(enz.columns):
        missing = required_enz_cols - set(enz.columns)
        raise KeyError(f"df_enzymes missing columns: {missing}")

    # Extract gene id (text after the last dot)
    pax["gene"] = (
        pax["string_external_id"]
        .astype(str)
        .str.split(".")
        .str[-1]
        .str.strip()
    )

    # Make abundance numeric and drop unusable rows
    pax["abundance"] = pd.to_numeric(pax["abundance"], errors="coerce")
    pax = pax.dropna(subset=["gene", "abundance"])

    # Aggregate in case there are multiple entries per gene (mean by default)
    pax_gene = (
        pax.groupby("gene", as_index=False)["abundance"]
        .mean()
        .rename(columns={"abundance": "protein_ppm"})
    )

    # Merge onto enzymes (left join to keep all enzymes)
    enz_mapped = enz.merge(pax_gene, on="gene", how="left")
    
    # Calculate molecular weights for each protein sequence # g/mol
    enz_mapped["molecular_weight"] = enz_mapped["sequence"].apply(calculate_molecular_weight)
    
    # Fraction of protein in the cell
    enz_mapped["protein_fraction"] = enz_mapped["protein_ppm"] / 1000000
    
    # Calculate protein_mol_gdcw: protein_fraction * p_total / molecular_weight
    # where p_total is the total protein content in g/gDCW
    # molecular_weight is in g/mol
    enz_mapped["protein_mol_gdcw"] = (
        enz_mapped["protein_fraction"] * p_total / enz_mapped["molecular_weight"]
    )
    
    # mol/gDCW to mmol/gDCW
    enz_mapped["protein_mmol_gdcw"] = enz_mapped["protein_mol_gdcw"] * 1000

    return enz_mapped


def process_enzyme_protein_mapping(enzyme_info_dfs: dict, paxdb_path: str, p_total: list):
    """
    Apply PaxDB protein mapping across all enzyme info dataframes and p_total values.
    
    Parameters:
        enzyme_info_dfs: dict
            Dictionary with condition names as keys and enzyme dataframes as values
        paxdb_path: str
            Path to PaxDB TSV file
        p_total: list
            List of total protein content values (g/gDCW) to test
    
    Returns:
        enzyme_protein_info_dfs: dict
            Nested dictionary structure:
            {condition_name: {p_total_value: mapped_dataframe}}
    """
    
    # Load PaxDB data
    print(f"Loading PaxDB data from: {paxdb_path}")
    paxdb_df = pd.read_csv(paxdb_path, sep="\t", comment="#", header=None, names=["gene_name", "string_external_id", "abundance"])
    print(f"PaxDB data loaded: {len(paxdb_df)} rows")
    
    # Initialize output dictionary
    enzyme_protein_info_dfs = {}
    
    # Get condition names and p_total values for progress tracking
    total_combinations = len(enzyme_info_dfs) * len(p_total)
    current_combination = 0
    
    print(f"Processing {len(enzyme_info_dfs)} conditions × {len(p_total)} p_total values = {total_combinations} combinations")
    
    # Double nested loop: for each condition, for each p_total value
    for condition_name, enzyme_df in enzyme_info_dfs.items():
        print(f"\nProcessing condition: {condition_name}")
        
        # Initialize nested dictionary for this condition
        enzyme_protein_info_dfs[condition_name] = {}
        
        for p_value in p_total:
            current_combination += 1
            print(f"  Processing p_total={p_value} ({current_combination}/{total_combinations})")
            
            # Apply map_paxdb_to_gene function
            try:
                mapped_df = map_paxdb_to_gene(
                    paxdb_df=paxdb_df,
                    df_enzymes=enzyme_df,
                    p_total=p_value
                )
                
                # Store the result
                enzyme_protein_info_dfs[condition_name][p_value] = mapped_df
                
                print(f"    Success: {len(mapped_df)} rows, {mapped_df['protein_ppm'].notna().sum()} with protein data")
                
            except Exception as e:
                print(f"    Error processing {condition_name} with p_total={p_value}: {str(e)}")
                # Store None or empty dataframe for failed combinations
                enzyme_protein_info_dfs[condition_name][p_value] = None
    
    print("\nCompleted processing all combinations")
    print(f"Results structure: {len(enzyme_protein_info_dfs)} conditions × {len(p_total)} p_total values")
    
    return enzyme_protein_info_dfs



def calculate_kapp_homomeric(enzyme_protein_info_dfs: dict):
    """
    Calculate kapp for homomeric enzymes for each condition and p_total combination.
    
    Parameters:
        enzyme_protein_info_dfs: dict
            Nested dictionary with structure {condition: {p_total: dataframe}}
    
    Returns:
        dict: Same structure as input but with added 'kcat_app' column in each dataframe
    """
    
    # Initialize output dictionary
    kapp_results = {}
    
    # Double nested loop: for each condition, for each p_total value
    for condition_name, p_total_dict in enzyme_protein_info_dfs.items():
        print(f"\nProcessing condition: {condition_name}")
        
        # Initialize nested dictionary for this condition
        kapp_results[condition_name] = {}
        
        for p_total_value, enzyme_df in p_total_dict.items():
            print(f"  Processing p_total={p_total_value}")
            
            # Skip if dataframe is None (failed processing)
            if enzyme_df is None:
                print("    Skipping - no data available")
                kapp_results[condition_name][p_total_value] = None
                continue
            
            # Work with a copy to avoid modifying original
            df_copy = enzyme_df.copy()
            
            print(f'    Rows before filtering homomeric: {len(df_copy)}')
            
            # Keep only homomeric enzymes
            df_copy = df_copy[df_copy['gpr_class'] == 'simple']
            print(f'    Rows after filtering homomeric: {len(df_copy)}')
            
            # Drop duplicate enzymes
            print(f'    Rows before filtering duplicates: {len(df_copy)}')
            df_copy = df_copy.drop_duplicates(subset=["gene", "SMILES"])
            print(f'    Rows after filtering duplicates: {len(df_copy)}')

            # Convert negative fluxes to positive (take absolute value)
            df_copy['flux_value'] = df_copy['flux_value'].abs()
            
            # COBRA fluxes are in mmol/gDW*h, convert to mmol/gDW*s
            df_copy['flux_value_per_sec'] = df_copy['flux_value'] / 3600  # mmol/gDW*s
            
            # Calculate kcat_app: flux (mmol/gDW*s) / enzyme concentration (mmol/gDCW) = kcat (1/s)
            # Handle division by zero by replacing with NaN
            df_copy['kcat_app'] = df_copy['flux_value_per_sec'] / df_copy['protein_mmol_gdcw']
            
            # Replace infinite values with NaN
            df_copy['kcat_app'] = df_copy['kcat_app'].replace([float('inf'), float('-inf')], float('nan'))
            
            # Count valid kcat_app values
            valid_kcat = df_copy['kcat_app'].notna().sum()
            print(f"    Calculated kcat_app for {valid_kcat} enzymes")
                
           
            # Store the processed dataframe
            kapp_results[condition_name][p_total_value] = df_copy
    
    print("\nCompleted kcat_app calculation for all conditions and p_total values")
    return kapp_results

def evaluate_kapp_homomeric(kapp_results: dict, upper_threshold: float = 1e6, lower_threshold: float = 1e-5):
    """
    Evaluate kapp for homomeric enzymes by filtering out unrealistic high and low values.
    
    Parameters:
        kapp_results: dict
            Nested dictionary with structure {condition: {p_total: dataframe}}
        upper_threshold: float
            Upper threshold for filtering kcat_app values (default: 1e6 s⁻¹)
            Rows with kcat_app > upper_threshold will be removed
        lower_threshold: float
            Lower threshold for filtering kcat_app values (default: 1e-4 s⁻¹)
            Rows with kcat_app < lower_threshold will be removed
    
    Returns:
        kapp_filtered_results: dict
            Same structure as input but with filtered dataframes
    """
    
    print(f"Filtering kcat_app values outside range: {lower_threshold:.0e} to {upper_threshold:.0e} s⁻¹")
    
    # Initialize output dictionary
    kapp_filtered_results = {}
    
    # Track filtering statistics
    total_original_rows = 0
    total_filtered_rows = 0
    total_removed_high = 0
    total_removed_low = 0
    
    # Double nested loop: for each condition, for each p_total value
    for condition_name, p_total_dict in kapp_results.items():
        print(f"\nProcessing condition: {condition_name}")
        
        # Initialize nested dictionary for this condition
        kapp_filtered_results[condition_name] = {}
        
        for p_total_value, df in p_total_dict.items():
            print(f"  Processing p_total={p_total_value}")
            
            # Skip if dataframe is None
            if df is None:
                print("    Skipping - no data available")
                kapp_filtered_results[condition_name][p_total_value] = None
                continue
            
            # Work with a copy to avoid modifying original
            df_filtered = df.copy()
            
            original_count = len(df_filtered)
            print(f"    Original rows: {original_count}")
            
            # Count values that will be removed for each threshold
            high_values = df_filtered['kcat_app'] > upper_threshold
            low_values = df_filtered['kcat_app'] < lower_threshold
            removed_high_count = high_values.sum()
            removed_low_count = low_values.sum()
            
            # Filter out rows where kcat_app is outside the acceptable range
            # Keep rows where kcat_app is NaN, or within [lower_threshold, upper_threshold]
            mask = (df_filtered['kcat_app'].isna()) | ((df_filtered['kcat_app'] >= lower_threshold) & (df_filtered['kcat_app'] <= upper_threshold))
            df_filtered = df_filtered[mask]
            
            filtered_count = len(df_filtered)
            total_removed_count = original_count - filtered_count
            
            print(f"    Filtered rows: {filtered_count}")
            print(f"    Removed total: {total_removed_count} (high: {removed_high_count}, low: {removed_low_count})")
            
            if removed_high_count > 0:
                high_removed_values = df[high_values]['kcat_app'].dropna()
                if len(high_removed_values) > 0:
                    print(f"    Removed high values range: {high_removed_values.min():.2e} to {high_removed_values.max():.2e} s⁻¹")
            
            if removed_low_count > 0:
                low_removed_values = df[low_values]['kcat_app'].dropna()
                if len(low_removed_values) > 0:
                    print(f"    Removed low values range: {low_removed_values.min():.2e} to {low_removed_values.max():.2e} s⁻¹")
            
            # Update statistics
            total_original_rows += original_count
            total_filtered_rows += filtered_count
            total_removed_high += removed_high_count
            total_removed_low += removed_low_count
            
            # Store the filtered dataframe
            kapp_filtered_results[condition_name][p_total_value] = df_filtered
    
    # Print summary statistics
    total_removed_rows = total_removed_high + total_removed_low
    print("\nFiltering Summary:")
    print(f"Total original rows: {total_original_rows}")
    print(f"Total filtered rows: {total_filtered_rows}")
    print(f"Total removed rows: {total_removed_rows}")
    print(f"  - Removed high values (>{upper_threshold:.0e}): {total_removed_high}")
    print(f"  - Removed low values (<{lower_threshold:.0e}): {total_removed_low}")
    if total_original_rows > 0:
        removal_percentage = (total_removed_rows / total_original_rows) * 100
        high_percentage = (total_removed_high / total_original_rows) * 100
        low_percentage = (total_removed_low / total_original_rows) * 100
        print(f"Percentage removed: {removal_percentage:.1f}% (high: {high_percentage:.1f}%, low: {low_percentage:.1f}%)")
    
    return kapp_filtered_results

def get_kmax_homomeric(kapp_results: dict):
    """
    Get the maximum kapp of each enzyme-substrate pair across all conditions and p_total combinations.
    
    Parameters:
        kapp_results: dict
            Nested dictionary with structure {condition: {p_total: dataframe}}
    Returns:
        kmax_results: pd.DataFrame
            DataFrame with columns: ['sequence', 'SMILES', 'kcat_app_max', 'condition_max', 'p_total_max']
            containing the maximum kcat_app value for each enzyme-substrate pair
    """
    
    print("Starting kmax analysis across all conditions and p_total values...")
    
    # List to collect all dataframes with metadata
    all_dataframes = []
    
    # Flatten the nested dictionary structure
    for condition_name, p_total_dict in kapp_results.items():
        for p_total_value, df in p_total_dict.items():
            # Skip None dataframes
            if df is None or len(df) == 0:
                continue
                
            # Add metadata columns to track source
            df_with_metadata = df.copy()
            df_with_metadata['source_condition'] = condition_name
            df_with_metadata['source_p_total'] = p_total_value
            
            # Only keep rows with valid kcat_app values
            df_with_metadata = df_with_metadata[df_with_metadata['kcat_app'].notna()]
            
            if len(df_with_metadata) > 0:
                all_dataframes.append(df_with_metadata)
                print(f"  Added {len(df_with_metadata)} valid entries from {condition_name}, p_total={p_total_value}")
    
    if not all_dataframes:
        print("No valid data found across all conditions")
        return pd.DataFrame(columns=['sequence', 'SMILES', 'kcat_app_max', 'condition_max', 'p_total_max'])
    
    # Concatenate all dataframes
    combined_df = pd.concat(all_dataframes, ignore_index=True)
    print(f"Combined dataframe has {len(combined_df)} total entries")
    
    # Group by enzyme-substrate pair (sequence + SMILES) and find maximum kcat_app
    print("Finding maximum kcat_app for each enzyme-substrate pair...")
    
    # Group by sequence and SMILES, then find the row with maximum kcat_app for each group
    kmax_results = (
        combined_df.loc[combined_df.groupby(['sequence', 'SMILES'])['kcat_app'].idxmax()]
        .reset_index(drop=True)
    )
    
    # Select and rename relevant columns for output
    output_columns = [
        'sequence', 'SMILES', 'kcat_app', 'source_condition', 'source_p_total',
        'gene', 'rxn', 'flux_value', 'FVA_upper', 'FVA_lower', 'protein_mmol_gdcw', 'subsystem' 
        # Additional useful columns
    ]
    
    # Keep only columns that exist in the dataframe
    available_columns = [col for col in output_columns if col in kmax_results.columns]
    kmax_results = kmax_results[available_columns].copy()
    
    # Rename columns for clarity
    column_renames = {
        'kcat_app': 'kcat_app_max',
        'source_condition': 'condition_max',
        'source_p_total': 'p_total_max'
    }
    kmax_results = kmax_results.rename(columns=column_renames)
    
    # Sort by kcat_app_max in descending order
    kmax_results = kmax_results.sort_values('kcat_app_max', ascending=False).reset_index(drop=True)
    
    print(f"Found maximum kcat_app values for {len(kmax_results)} unique enzyme-substrate pairs")
    print(f"kcat_app_max range: {kmax_results['kcat_app_max'].min():.2e} to {kmax_results['kcat_app_max'].max():.2e} s⁻¹")
    
    # Show summary statistics
    condition_counts = kmax_results['condition_max'].value_counts()
    print("Maximum values found across conditions:")
    for condition, count in condition_counts.items():
        print(f"  {condition}: {count} enzyme-substrate pairs")
    
    return kmax_results


def get_eta(kapp_results: dict, kmax_results: pd.DataFrame):
    """
    Calculate eta (kapp/kmax) for each enzyme-substrate pair across all conditions and p_total values.
    
    Parameters:
        kapp_results: dict
            Nested dictionary with structure {condition: {p_total: dataframe}}
        kmax_results: pd.DataFrame
            DataFrame with maximum kcat_app values for each enzyme-substrate pair
    
    Returns:
        tuple: (kapp_results_with_eta, kmax_results_with_variance)
            - kapp_results_with_eta: dict with same structure as input but with 'eta' column added
            - kmax_results_with_variance: DataFrame with added variance columns (eta_mean, eta_stdev, eta_min, eta_max, eta_cv)
    """
    
    print("Calculating eta (kapp/kmax) for all conditions...")
    
    # Initialize output dictionary
    kapp_results_with_eta = {}
    
    # List to collect all eta values for variance calculation
    all_eta_values = []
    
    # Process each condition and p_total combination
    for condition_name, p_total_dict in kapp_results.items():
        print(f"\nProcessing condition: {condition_name}")
        
        # Initialize nested dictionary for this condition
        kapp_results_with_eta[condition_name] = {}
        
        for p_total_value, df in p_total_dict.items():
            print(f"  Processing p_total={p_total_value}")
            
            # Skip if dataframe is None
            if df is None:
                print("    Skipping - no data available")
                kapp_results_with_eta[condition_name][p_total_value] = None
                continue
            
            # Work with a copy
            df_with_eta = df.copy()
            
            # Merge with kmax_results to get the maximum kcat_app value
            df_with_eta = pd.merge(
                df_with_eta,
                kmax_results[['sequence', 'SMILES', 'kcat_app_max']],
                on=['sequence', 'SMILES'],
                how='left'
            )
            
            # Calculate eta = kcat_app / kcat_app_max
            df_with_eta['eta'] = df_with_eta['kcat_app'] / df_with_eta['kcat_app_max']
            
            # Replace infinite values with NaN
            df_with_eta['eta'] = df_with_eta['eta'].replace([float('inf'), float('-inf')], float('nan'))
            
            # Count valid eta values
            valid_eta = df_with_eta['eta'].notna().sum()
            print(f"    Calculated eta for {valid_eta} enzyme-substrate pairs")
            
            # Store the dataframe with eta
            kapp_results_with_eta[condition_name][p_total_value] = df_with_eta
            
            # Collect eta values for variance calculation
            eta_data = df_with_eta[['sequence', 'SMILES', 'eta']].copy()
            eta_data['source_condition'] = condition_name
            eta_data['source_p_total'] = p_total_value
            eta_data = eta_data[eta_data['eta'].notna()]  # Keep only valid eta values
            
            if len(eta_data) > 0:
                all_eta_values.append(eta_data)
    
    # Calculate variance metrics for each enzyme-substrate pair
    print("\nCalculating variance metrics for each enzyme-substrate pair...")
    
    if not all_eta_values:
        print("No valid eta values found")
        # Add empty variance columns to kmax_results
        kmax_with_variance = kmax_results.copy()
        kmax_with_variance['eta_mean'] = float('nan')
        kmax_with_variance['eta_stdev'] = float('nan')
        kmax_with_variance['eta_min'] = float('nan')
        kmax_with_variance['eta_max'] = float('nan')
        kmax_with_variance['eta_cv'] = float('nan')
        return kapp_results_with_eta, kmax_with_variance
    
    # Concatenate all eta values
    all_eta_df = pd.concat(all_eta_values, ignore_index=True)
    print(f"Collected {len(all_eta_df)} eta values across all conditions")
    
    # Group by enzyme-substrate pair and calculate variance metrics
    variance_metrics = all_eta_df.groupby(['sequence', 'SMILES'])['eta'].agg([
        ('eta_mean', 'mean'),
        ('eta_stdev', 'std'),
        ('eta_min', 'min'),
        ('eta_max', 'max'),
        ('eta_count', 'count')
    ]).reset_index()
    
    # Calculate coefficient of variation (CV = stdev / mean)
    variance_metrics['eta_cv'] = variance_metrics['eta_stdev'] / variance_metrics['eta_mean']
    
    # Replace infinite CV values with NaN
    variance_metrics['eta_cv'] = variance_metrics['eta_cv'].replace([float('inf'), float('-inf')], float('nan'))
    
    # Merge variance metrics with kmax_results
    kmax_with_variance = pd.merge(
        kmax_results,
        variance_metrics[['sequence', 'SMILES', 'eta_mean', 'eta_stdev', 'eta_min', 'eta_max', 'eta_cv']],
        on=['sequence', 'SMILES'],
        how='left'
    )
    
    print(f"Added variance metrics to {len(kmax_with_variance)} enzyme-substrate pairs")
    print(f"Pairs with variance data: {kmax_with_variance['eta_mean'].notna().sum()}")
    
    # Show summary statistics
    if len(variance_metrics) > 0:
        print("\nEta variance summary:")
        print(f"  Mean eta_mean: {variance_metrics['eta_mean'].mean():.3f}")
        print(f"  Mean eta_stdev: {variance_metrics['eta_stdev'].mean():.3f}")
        print(f"  Mean eta_cv: {variance_metrics['eta_cv'].mean():.3f}")
    
    print("\nCompleted eta calculation for all conditions")
    return kapp_results_with_eta, kmax_with_variance


def identify_heteromeric_enzymes(enzyme_info_dfs: dict):
    """
    Filter and prepare heteromeric enzyme data from the enzyme_info_dfs.
    
    Parameters:
        enzyme_info_dfs: dict
            Dictionary with condition names as keys and enzyme dataframes as values
    
    Returns:
        tuple: (heteromeric_enzyme_dfs, heteromeric_metadata)
            - heteromeric_enzyme_dfs: dict with same structure containing only heteromeric enzymes
            - heteromeric_metadata: dict with information about each complex
    """
    
    print("Identifying heteromeric enzyme complexes...")
    
    # Initialize output dictionaries
    heteromeric_enzyme_dfs = {}
    heteromeric_metadata = {}
    
    # Process each condition
    for condition_name, df in enzyme_info_dfs.items():
        print(f"\nProcessing condition: {condition_name}")
        
        # Filter for heteromeric complexes (gpr_class == 'and_only')
        # Only keep heteromeric complexes that are connected by just AND logic
        # OR and AND cases will be handled separately
        heteromeric_df = df[df['gpr_class'] == 'and_only'].copy()
        
        print(f"  Found {len(heteromeric_df)} heteromeric enzyme entries")
        
        if len(heteromeric_df) == 0:
            heteromeric_enzyme_dfs[condition_name] = heteromeric_df
            continue
        
        # Group by reaction to identify complexes
        for rxn_id, group in heteromeric_df.groupby('rxn'):
            complex_id = f"{rxn_id}"
            
            # Get all genes (subunits) for this complex
            genes = group['gene'].unique().tolist()
            
            # Check if all subunits have sequence data
            missing_sequences = group['sequence'].isna().sum()
            all_sequences_available = (missing_sequences == 0)
            
            # Store metadata
            if complex_id not in heteromeric_metadata:
                heteromeric_metadata[complex_id] = {
                    'complex_id': complex_id,
                    'reaction_id': rxn_id,
                    'genes': genes,
                    'n_subunits': len(genes),
                    'all_sequences_available': all_sequences_available,
                    'missing_sequences_count': missing_sequences
                }
        
        # Store the heteromeric dataframe for this condition
        heteromeric_enzyme_dfs[condition_name] = heteromeric_df
        
        print(f"  Identified {len(heteromeric_df['rxn'].unique())} unique heteromeric complexes")
        print(f"  Total subunit entries: {len(heteromeric_df)}")
    
    # Print summary of metadata
    print("\n" + "="*60)
    print("HETEROMERIC COMPLEXES SUMMARY")
    print("="*60)
    print(f"Total unique complexes identified: {len(heteromeric_metadata)}")
    
    complexes_with_all_sequences = sum(1 for v in heteromeric_metadata.values() if v['all_sequences_available'])
    print(f"Complexes with all sequences available: {complexes_with_all_sequences}")
    print(f"Complexes with missing sequences: {len(heteromeric_metadata) - complexes_with_all_sequences}")
    
    return heteromeric_enzyme_dfs, heteromeric_metadata


def map_heteromeric_subunits_to_proteomics(heteromeric_enzyme_dfs: dict, 
                                           paxdb_path: str, 
                                           p_total: list):
    """
    Map all subunits of heteromeric complexes to proteomics data and aggregate.
    
    Parameters:
        heteromeric_enzyme_dfs: dict
            Dictionary with condition names as keys and heteromeric enzyme dataframes as values
        paxdb_path: str
            Path to PaxDB TSV file
        p_total: list
            List of total protein content values (g/gDCW) to test
    
    Returns:
        heteromeric_protein_info_dfs: dict
            Nested dictionary {condition: {p_total: dataframe}}
            Each row represents a complete enzyme complex with aggregated subunit data
    """
    
    # Load PaxDB data
    print(f"Loading PaxDB data from: {paxdb_path}")
    paxdb_df = pd.read_csv(paxdb_path, sep="\t", comment="#", header=None, 
                           names=["gene_name", "string_external_id", "abundance"])
    
    # Extract gene id from PaxDB
    paxdb_df["gene"] = (
        paxdb_df["string_external_id"]
        .astype(str)
        .str.split(".")
        .str[-1]
        .str.strip()
    )
    paxdb_df["abundance"] = pd.to_numeric(paxdb_df["abundance"], errors="coerce")
    paxdb_df = paxdb_df.dropna(subset=["gene", "abundance"])
    
    # Aggregate by gene
    paxdb_gene = (
        paxdb_df.groupby("gene", as_index=False)["abundance"]
        .mean()
        .rename(columns={"abundance": "protein_ppm"})
    )
    
    print(f"PaxDB data loaded: {len(paxdb_gene)} unique genes with abundance data")
    
    # Initialize output dictionary
    heteromeric_protein_info_dfs = {}
    
    total_combinations = len(heteromeric_enzyme_dfs) * len(p_total)
    current_combination = 0
    
    print(f"\nProcessing {len(heteromeric_enzyme_dfs)} conditions × {len(p_total)} p_total values")
    
    # Process each condition and p_total combination
    for condition_name, heteromeric_df in heteromeric_enzyme_dfs.items():
        print(f"\nProcessing condition: {condition_name}")
        
        heteromeric_protein_info_dfs[condition_name] = {}
        
        for p_value in p_total:
            current_combination += 1
            print(f"  Processing p_total={p_value} ({current_combination}/{total_combinations})")
            
            if heteromeric_df is None or len(heteromeric_df) == 0:
                print("    No heteromeric enzymes in this condition")
                heteromeric_protein_info_dfs[condition_name][p_value] = None
                continue
            
            # Work with a copy
            df_copy = heteromeric_df.copy()
            
            # Calculate molecular weight for each subunit
            df_copy['molecular_weight'] = df_copy['sequence'].apply(calculate_molecular_weight)
            
            # Merge with PaxDB abundance data
            df_copy = pd.merge(df_copy, paxdb_gene, on='gene', how='left')
            
            # Calculate protein fraction and concentration for each subunit
            df_copy['protein_fraction'] = df_copy['protein_ppm'] / 1e6
            df_copy['protein_mol_gdcw'] = (
                df_copy['protein_fraction'] * p_value / df_copy['molecular_weight']
            )
            df_copy['protein_mmol_gdcw'] = df_copy['protein_mol_gdcw'] * 1000
            
            # Calculate mass per subunit (g/gDCW) = protein_fraction * p_total
            df_copy['subunit_g_gdcw'] = df_copy['protein_fraction'] * p_value
            
            # Group by reaction to aggregate subunits into complexes
            complex_rows = []
            
            for rxn_id, group in df_copy.groupby('rxn'):
                # Check if all subunits have abundance data
                missing_abundance = group['protein_ppm'].isna().sum()
                
                if missing_abundance > 0:
                    # Skip complexes with missing subunit data
                    continue
                
                # Aggregate subunit information
                complex_genes = sorted(group['gene'].unique().tolist())
                complex_sequences = group['sequence'].tolist()
                
                # Sum molecular weights of all subunits
                complex_mw = group['molecular_weight'].sum()
                
                # Sum mass fractions of all subunits (Davidi et al. Eq. S3)
                complex_g_gdcw = group['subunit_g_gdcw'].sum()
                
                # Create a dictionary with subunit details
                subunit_details = []
                for _, subunit in group.iterrows():
                    subunit_details.append({
                        'gene': subunit['gene'],
                        'protein_ppm': subunit['protein_ppm'],
                        'molecular_weight': subunit['molecular_weight'],
                        'g_gdcw': subunit['subunit_g_gdcw']
                    })
                
                # Create complex row (take first row as template and modify)
                complex_row = group.iloc[0].copy()
                complex_row['complex_genes'] = ','.join(complex_genes)
                complex_row['n_subunits'] = len(complex_genes)
                complex_row['complex_molecular_weight'] = complex_mw
                complex_row['complex_g_gdcw'] = complex_g_gdcw
                complex_row['subunit_details'] = json.dumps(subunit_details)
                
                # Store gene list as a string for identification
                complex_row['gene'] = ','.join(complex_genes)
                
                complex_rows.append(complex_row)
            
            if len(complex_rows) == 0:
                print("    No complete complexes with all subunit data")
                heteromeric_protein_info_dfs[condition_name][p_value] = pd.DataFrame()
            else:
                result_df = pd.DataFrame(complex_rows)
                heteromeric_protein_info_dfs[condition_name][p_value] = result_df
                print(f"    Successfully processed {len(result_df)} complete complexes")
    
    print("\nCompleted heteromeric protein mapping")
    return heteromeric_protein_info_dfs


def calculate_specific_activity_heteromeric(heteromeric_protein_info_dfs: dict):
    """
    Calculate specific activity (SA_app) for heteromeric complexes.
    
    Parameters:
        heteromeric_protein_info_dfs: dict
            Nested dictionary {condition: {p_total: dataframe}}
    
    Returns:
        dict: Same structure with added 'SA_app' column (μmol/mg/min)
    """
    
    print("Calculating specific activity for heteromeric complexes...")
    
    # Initialize output dictionary
    SA_results = {}
    
    # Process each condition and p_total combination
    for condition_name, p_total_dict in heteromeric_protein_info_dfs.items():
        print(f"\nProcessing condition: {condition_name}")
        
        SA_results[condition_name] = {}
        
        for p_total_value, df in p_total_dict.items():
            print(f"  Processing p_total={p_total_value}")
            
            # Skip if dataframe is None or empty
            if df is None or len(df) == 0:
                print("    Skipping - no data available")
                SA_results[condition_name][p_total_value] = None
                continue
            
            # Work with a copy
            df_copy = df.copy()
            
            print(f"    Processing {len(df_copy)} heteromeric complexes")
            
            # Drop duplicate complexes (same gene set and substrate)
            df_copy = df_copy.drop_duplicates(subset=['complex_genes', 'SMILES'])
            print(f"    After removing duplicates: {len(df_copy)} complexes")
            
            # Convert flux to absolute value
            df_copy['flux_value'] = df_copy['flux_value'].abs()
            
            # Convert flux from mmol/gDCW/h to μmol/gDCW/min
            # mmol/gDCW/h * 1000 (to μmol) / 60 (to min) = μmol/gDCW/min
            df_copy['flux_umol_gdcw_min'] = df_copy['flux_value'] * 1000 / 60
            
            # Calculate specific activity: flux (μmol/gDCW/min) / complex concentration (mg/gDCW)
            # SA_app units: μmol/mg/min
            # complex_g_gdcw is in g/gDCW, convert to mg/gDCW by multiplying by 1000
            df_copy['SA_app'] = df_copy['flux_umol_gdcw_min'] / (df_copy['complex_g_gdcw'] * 1000)
            
            # Replace infinite values with NaN
            df_copy['SA_app'] = df_copy['SA_app'].replace([float('inf'), float('-inf')], float('nan'))
            
            # Count valid SA_app values
            valid_SA = df_copy['SA_app'].notna().sum()
            print(f"    Calculated SA_app for {valid_SA} complexes")
            
            # Store the processed dataframe
            SA_results[condition_name][p_total_value] = df_copy
    
    print("\nCompleted SA_app calculation for all conditions and p_total values")
    return SA_results


def evaluate_specific_activity_heteromeric(SA_results: dict, 
                                           upper_threshold: float = 1e4, 
                                           lower_threshold: float = 1e-3):
    """
    Filter unrealistic specific activity values for heteromeric complexes.
    
    Parameters:
        SA_results: dict
            Nested dictionary {condition: {p_total: dataframe}}
        upper_threshold: float
            Upper threshold for SA_app (default: 1e4 μmol/mg/min)
        lower_threshold: float
            Lower threshold for SA_app (default: 1e-3 μmol/mg/min)
    
    Returns:
        SA_filtered_results: dict
            Same structure with filtered dataframes
    """
    
    print(f"Filtering SA_app values outside range: {lower_threshold:.0e} to {upper_threshold:.0e} μmol/mg/min")
    
    # Initialize output dictionary
    SA_filtered_results = {}
    
    # Track filtering statistics
    total_original_rows = 0
    total_filtered_rows = 0
    total_removed_high = 0
    total_removed_low = 0
    
    # Process each condition and p_total
    for condition_name, p_total_dict in SA_results.items():
        print(f"\nProcessing condition: {condition_name}")
        
        SA_filtered_results[condition_name] = {}
        
        for p_total_value, df in p_total_dict.items():
            print(f"  Processing p_total={p_total_value}")
            
            # Skip if dataframe is None
            if df is None or len(df) == 0:
                print("    Skipping - no data available")
                SA_filtered_results[condition_name][p_total_value] = None
                continue
            
            # Work with a copy
            df_filtered = df.copy()
            
            original_count = len(df_filtered)
            print(f"    Original rows: {original_count}")
            
            # Count values that will be removed
            high_values = df_filtered['SA_app'] > upper_threshold
            low_values = df_filtered['SA_app'] < lower_threshold
            removed_high_count = high_values.sum()
            removed_low_count = low_values.sum()
            
            # Filter
            mask = (df_filtered['SA_app'].isna()) | \
                   ((df_filtered['SA_app'] >= lower_threshold) & (df_filtered['SA_app'] <= upper_threshold))
            df_filtered = df_filtered[mask]
            
            filtered_count = len(df_filtered)
            total_removed_count = original_count - filtered_count
            
            print(f"    Filtered rows: {filtered_count}")
            print(f"    Removed total: {total_removed_count} (high: {removed_high_count}, low: {removed_low_count})")
            
            # Update statistics
            total_original_rows += original_count
            total_filtered_rows += filtered_count
            total_removed_high += removed_high_count
            total_removed_low += removed_low_count
            
            # Store the filtered dataframe
            SA_filtered_results[condition_name][p_total_value] = df_filtered
    
    # Print summary
    total_removed_rows = total_removed_high + total_removed_low
    print("\nFiltering Summary:")
    print(f"Total original rows: {total_original_rows}")
    print(f"Total filtered rows: {total_filtered_rows}")
    print(f"Total removed rows: {total_removed_rows}")
    print(f"  - Removed high values (>{upper_threshold:.0e}): {total_removed_high}")
    print(f"  - Removed low values (<{lower_threshold:.0e}): {total_removed_low}")
    if total_original_rows > 0:
        removal_percentage = (total_removed_rows / total_original_rows) * 100
        print(f"Percentage removed: {removal_percentage:.1f}%")
    
    return SA_filtered_results


def get_SA_max_heteromeric(SA_results: dict):
    """
    Find maximum specific activity across all conditions for heteromeric complexes.
    
    Parameters:
        SA_results: dict
            Nested dictionary {condition: {p_total: dataframe}}
    
    Returns:
        pd.DataFrame: DataFrame with SA_max for each complex-substrate pair
    """
    
    print("Finding SA_max for heteromeric complexes...")
    
    # List to collect all dataframes
    all_dataframes = []
    
    # Flatten nested structure
    for condition_name, p_total_dict in SA_results.items():
        for p_total_value, df in p_total_dict.items():
            if df is None or len(df) == 0:
                continue
            
            # Add metadata
            df_with_metadata = df.copy()
            df_with_metadata['source_condition'] = condition_name
            df_with_metadata['source_p_total'] = p_total_value
            
            # Keep only valid SA_app values
            df_with_metadata = df_with_metadata[df_with_metadata['SA_app'].notna()]
            
            if len(df_with_metadata) > 0:
                all_dataframes.append(df_with_metadata)
                print(f"  Added {len(df_with_metadata)} entries from {condition_name}, p_total={p_total_value}")
    
    if not all_dataframes:
        print("No valid data found")
        return pd.DataFrame(columns=['complex_genes', 'SMILES', 'SA_app_max', 'condition_max', 'p_total_max'])
    
    # Concatenate
    combined_df = pd.concat(all_dataframes, ignore_index=True)
    print(f"Combined dataframe has {len(combined_df)} total entries")
    
    # Find maximum SA_app for each complex-substrate pair
    print("Finding maximum SA_app for each complex-substrate pair...")
    
    SA_max_results = (
        combined_df.loc[combined_df.groupby(['complex_genes', 'SMILES'])['SA_app'].idxmax()]
        .reset_index(drop=True)
    )
    
    # Select relevant columns
    output_columns = [
        'complex_genes', 'n_subunits', 'SMILES', 'SA_app', 'source_condition', 'source_p_total',
        'rxn', 'flux_value', 'FVA_upper', 'FVA_lower', 'complex_g_gdcw', 'complex_molecular_weight',
        'subsystem', 'subunit_details'
    ]
    
    available_columns = [col for col in output_columns if col in SA_max_results.columns]
    SA_max_results = SA_max_results[available_columns].copy()
    
    # Rename columns
    column_renames = {
        'SA_app': 'SA_app_max',
        'source_condition': 'condition_max',
        'source_p_total': 'p_total_max'
    }
    SA_max_results = SA_max_results.rename(columns=column_renames)
    
    # Sort by SA_app_max
    SA_max_results = SA_max_results.sort_values('SA_app_max', ascending=False).reset_index(drop=True)
    
    print(f"Found maximum SA_app for {len(SA_max_results)} unique complex-substrate pairs")
    print(f"SA_app_max range: {SA_max_results['SA_app_max'].min():.2e} to {SA_max_results['SA_app_max'].max():.2e} μmol/mg/min")
    
    # Show summary
    condition_counts = SA_max_results['condition_max'].value_counts()
    print("Maximum values found across conditions:")
    for condition, count in condition_counts.items():
        print(f"  {condition}: {count} complex-substrate pairs")
    
    return SA_max_results


def get_eta_heteromeric(SA_results: dict, SA_max_results: pd.DataFrame):
    """
    Calculate eta (SA_app/SA_max) and variance metrics for heteromeric complexes.
    
    Parameters:
        SA_results: dict
            Nested dictionary {condition: {p_total: dataframe}}
        SA_max_results: pd.DataFrame
            DataFrame with SA_max values
    
    Returns:
        tuple: (SA_results_with_eta, SA_max_results_with_variance)
    """
    
    print("Calculating eta (SA_app/SA_max) for heteromeric complexes...")
    
    # Initialize output dictionary
    SA_results_with_eta = {}
    
    # List to collect eta values
    all_eta_values = []
    
    # Process each condition and p_total
    for condition_name, p_total_dict in SA_results.items():
        print(f"\nProcessing condition: {condition_name}")
        
        SA_results_with_eta[condition_name] = {}
        
        for p_total_value, df in p_total_dict.items():
            print(f"  Processing p_total={p_total_value}")
            
            # Skip if dataframe is None
            if df is None or len(df) == 0:
                print("    Skipping - no data available")
                SA_results_with_eta[condition_name][p_total_value] = None
                continue
            
            # Work with a copy
            df_with_eta = df.copy()
            
            # Merge with SA_max_results
            df_with_eta = pd.merge(
                df_with_eta,
                SA_max_results[['complex_genes', 'SMILES', 'SA_app_max']],
                on=['complex_genes', 'SMILES'],
                how='left'
            )
            
            # Calculate eta
            df_with_eta['eta'] = df_with_eta['SA_app'] / df_with_eta['SA_app_max']
            
            # Replace infinite values with NaN
            df_with_eta['eta'] = df_with_eta['eta'].replace([float('inf'), float('-inf')], float('nan'))
            
            # Count valid eta values
            valid_eta = df_with_eta['eta'].notna().sum()
            print(f"    Calculated eta for {valid_eta} complex-substrate pairs")
            
            # Store
            SA_results_with_eta[condition_name][p_total_value] = df_with_eta
            
            # Collect eta values for variance
            eta_data = df_with_eta[['complex_genes', 'SMILES', 'eta']].copy()
            eta_data['source_condition'] = condition_name
            eta_data['source_p_total'] = p_total_value
            eta_data = eta_data[eta_data['eta'].notna()]
            
            if len(eta_data) > 0:
                all_eta_values.append(eta_data)
    
    # Calculate variance metrics
    print("\nCalculating variance metrics...")
    
    if not all_eta_values:
        print("No valid eta values found")
        SA_max_with_variance = SA_max_results.copy()
        SA_max_with_variance['eta_mean'] = float('nan')
        SA_max_with_variance['eta_stdev'] = float('nan')
        SA_max_with_variance['eta_min'] = float('nan')
        SA_max_with_variance['eta_max'] = float('nan')
        SA_max_with_variance['eta_cv'] = float('nan')
        return SA_results_with_eta, SA_max_with_variance
    
    # Concatenate
    all_eta_df = pd.concat(all_eta_values, ignore_index=True)
    print(f"Collected {len(all_eta_df)} eta values")
    
    # Calculate variance metrics
    variance_metrics = all_eta_df.groupby(['complex_genes', 'SMILES'])['eta'].agg([
        ('eta_mean', 'mean'),
        ('eta_stdev', 'std'),
        ('eta_min', 'min'),
        ('eta_max', 'max'),
        ('eta_count', 'count')
    ]).reset_index()
    
    # Calculate CV
    variance_metrics['eta_cv'] = variance_metrics['eta_stdev'] / variance_metrics['eta_mean']
    variance_metrics['eta_cv'] = variance_metrics['eta_cv'].replace([float('inf'), float('-inf')], float('nan'))
    
    # Merge with SA_max_results
    SA_max_with_variance = pd.merge(
        SA_max_results,
        variance_metrics[['complex_genes', 'SMILES', 'eta_mean', 'eta_stdev', 'eta_min', 'eta_max', 'eta_cv']],
        on=['complex_genes', 'SMILES'],
        how='left'
    )
    
    print(f"Added variance metrics to {len(SA_max_with_variance)} complex-substrate pairs")
    
    return SA_results_with_eta, SA_max_with_variance


def combine_homomeric_heteromeric_results(kmax_results: pd.DataFrame, 
                                          SA_max_results: pd.DataFrame):
    """
    Combine homomeric and heteromeric results into a single dataframe.
    
    Parameters:
        kmax_results: pd.DataFrame
            Homomeric kmax results
        SA_max_results: pd.DataFrame
            Heteromeric SA_max results
    
    Returns:
        pd.DataFrame: Combined results with enzyme_type column
    """
    
    print("Combining homomeric and heteromeric results...")
    
    # Add enzyme type column
    kmax_with_type = kmax_results.copy()
    kmax_with_type['enzyme_type'] = 'homomeric'
    kmax_with_type['enzyme_id'] = kmax_with_type['gene']
    kmax_with_type['rate_value'] = kmax_with_type['kcat_app_max']
    kmax_with_type['rate_units'] = 's^-1'
    
    SA_max_with_type = SA_max_results.copy()
    SA_max_with_type['enzyme_type'] = 'heteromeric'
    SA_max_with_type['enzyme_id'] = SA_max_with_type['complex_genes']
    SA_max_with_type['rate_value'] = SA_max_with_type['SA_app_max']
    SA_max_with_type['rate_units'] = 'μmol/mg/min'
    
    # Convert heteromeric SA to kcat if desired (optional)
    # kcat [s^-1] = SA [μmol/mg/min] * MW [mg/μmol] / 60 [s/min]
    # kcat = SA * MW / 60, where MW is in g/mol = mg/mmol = mg/(1000*μmol)
    if 'complex_molecular_weight' in SA_max_with_type.columns:
        SA_max_with_type['kcat_app_max_converted'] = (
            SA_max_with_type['SA_app_max'] * 
            SA_max_with_type['complex_molecular_weight'] / 
            1000 /  # Convert g/mol to mg/μmol
            60      # Convert minutes to seconds
        )
    
    # Select common columns for merging
    common_columns = [
        'enzyme_type', 'enzyme_id', 'rate_value', 'rate_units',
        'SMILES', 'rxn', 'condition_max', 'p_total_max',
        'flux_value', 'subsystem',
        'eta_mean', 'eta_stdev', 'eta_min', 'eta_max', 'eta_cv'
    ]
    
    # Only keep columns that exist in both dataframes
    kmax_available = [col for col in common_columns if col in kmax_with_type.columns]
    SA_available = [col for col in common_columns if col in SA_max_with_type.columns]
    
    # Add specific columns for each type
    kmax_specific = ['gene', 'sequence', 'kcat_app_max', 'protein_mmol_gdcw']
    SA_specific = ['complex_genes', 'n_subunits', 'SA_app_max', 'complex_g_gdcw', 'complex_molecular_weight']
    
    if 'kcat_app_max_converted' in SA_max_with_type.columns:
        SA_specific.append('kcat_app_max_converted')
    
    kmax_cols = kmax_available + [col for col in kmax_specific if col in kmax_with_type.columns]
    SA_cols = SA_available + [col for col in SA_specific if col in SA_max_with_type.columns]
    
    # Create subset dataframes
    kmax_subset = kmax_with_type[kmax_cols].copy()
    SA_subset = SA_max_with_type[SA_cols].copy()
    
    # Concatenate
    combined_results = pd.concat([kmax_subset, SA_subset], ignore_index=True, sort=False)
    
    # Sort by rate_value (descending)
    combined_results = combined_results.sort_values('rate_value', ascending=False).reset_index(drop=True)
    
    print(f"\nCombined results:")
    print(f"  Homomeric enzymes: {len(kmax_results)}")
    print(f"  Heteromeric complexes: {len(SA_max_results)}")
    print(f"  Total: {len(combined_results)}")
    
    # Summary statistics
    homomeric_count = (combined_results['enzyme_type'] == 'homomeric').sum()
    heteromeric_count = (combined_results['enzyme_type'] == 'heteromeric').sum()
    
    print(f"\nBreakdown by enzyme type:")
    print(f"  Homomeric: {homomeric_count} ({homomeric_count/len(combined_results)*100:.1f}%)")
    print(f"  Heteromeric: {heteromeric_count} ({heteromeric_count/len(combined_results)*100:.1f}%)")
    
    return combined_results


def expand_heteromeric_to_subunits(SA_max_results: pd.DataFrame):
    """
    Expand heteromeric complex results to create one row per subunit.
    Each subunit row inherits the SA_app_max from its parent complex.
    
    Parameters:
        SA_max_results: pd.DataFrame
            Heteromeric SA_max results with one row per complex
    
    Returns:
        pd.DataFrame: Expanded dataframe with columns:
            - sequence: Individual subunit sequence
            - gene: Individual subunit gene
            - SMILES: Substrate
            - SA_app_max: Specific activity (same for all subunits in a complex)
            - complex_genes: Original complex identifier
            - n_subunits: Number of subunits in complex
            - subunit_role: Position/identifier within complex
    """
    
    print("Expanding heteromeric complexes to individual subunits...")
    
    if SA_max_results is None or len(SA_max_results) == 0:
        print("  No heteromeric data to expand")
        return pd.DataFrame(columns=['sequence', 'gene', 'SMILES', 'SA_app_max'])
    
    expanded_rows = []
    
    for idx, row in SA_max_results.iterrows():
        # Parse subunit details from JSON
        try:
            subunit_details = json.loads(row['subunit_details'])
        except (KeyError, json.JSONDecodeError, TypeError):
            print(f"  Warning: Could not parse subunit_details for complex {row.get('complex_genes', 'unknown')}")
            continue
        
        # Get complex-level information
        complex_genes_list = row['complex_genes'].split(',')
        
        # For each subunit in the complex
        for i, subunit in enumerate(subunit_details):
            # Create a new row for this subunit
            subunit_row = row.copy()
            
            # Override with subunit-specific information
            subunit_row['gene'] = subunit['gene']
            subunit_row['subunit_role'] = f"subunit_{i+1}_of_{len(subunit_details)}"
            subunit_row['subunit_molecular_weight'] = subunit['molecular_weight']
            subunit_row['subunit_protein_ppm'] = subunit['protein_ppm']
            subunit_row['subunit_g_gdcw'] = subunit['g_gdcw']
            
            # Note: We need to get the sequence for this subunit
            # The sequence is not stored in subunit_details, so we'll need to handle this
            # We'll add a placeholder that should be filled by merging with sequence data
            
            expanded_rows.append(subunit_row)
    
    if not expanded_rows:
        print("  No subunits could be expanded")
        return pd.DataFrame(columns=['sequence', 'gene', 'SMILES', 'SA_app_max'])
    
    # Create expanded dataframe
    expanded_df = pd.DataFrame(expanded_rows)
    
    print(f"  Expanded {len(SA_max_results)} complexes into {len(expanded_df)} subunit rows")
    print(f"  Average subunits per complex: {len(expanded_df)/len(SA_max_results):.1f}")
    
    return expanded_df


def create_unified_kcat_SA(kmax_results: pd.DataFrame, 
                                           SA_max_results: pd.DataFrame,
                                           sequence_df: pd.DataFrame):
    """
    Create a unified dataframe with sequence-substrate pairs from both homomeric 
    and heteromeric enzymes, with appropriate rate metrics.
    
    For homomeric: sequence maps to kcat_app_max (s⁻¹)
    For heteromeric: each subunit sequence maps to SA_app_max (μmol/mg/min) of its complex
    
    Parameters:
        kmax_results: pd.DataFrame
            Homomeric kmax results
        SA_max_results: pd.DataFrame
            Heteromeric SA_max results
        sequence_df: str or pd.DataFrame
            Sequence dataframe or path to CSV file
    
    Returns:
        pd.DataFrame: Unified view with columns:
            - sequence: Protein sequence
            - gene: Gene identifier
            - SMILES: Substrate
            - rate_value: Either kcat_app_max or SA_app_max
            - rate_units: Either 's^-1' or 'μmol/mg/min'
            - enzyme_type: 'homomeric' or 'heteromeric'
            - molecular_weight: Protein MW (g/mol)
            - [for heteromeric: complex_genes, n_subunits, subunit_role]
    """
    
    print("Creating unified sequence-substrate view...")
    
    sequence_df = load_dataframe_if_path(sequence_df)
    
    # Process homomeric data
    homo_view = kmax_results.copy()
    homo_view['enzyme_type'] = 'homomeric'
    homo_view['rate_value'] = homo_view['kcat_app_max']
    homo_view['rate_units'] = 's^-1'
    
    # Select relevant columns for homomeric
    homo_columns = [
        'sequence', 'gene', 'SMILES', 'rate_value', 'rate_units', 'enzyme_type',
        'rxn', 'condition_max', 'p_total_max', 'flux_value', 'subsystem',
        'eta_mean', 'eta_stdev', 'eta_min', 'eta_max', 'eta_cv',
        'kcat_app_max', 'protein_mmol_gdcw'
    ]
    homo_columns = [col for col in homo_columns if col in homo_view.columns]
    homo_view = homo_view[homo_columns].copy()
    
    # Process heteromeric data - expand to subunits
    hetero_expanded = expand_heteromeric_to_subunits(SA_max_results)
    
    if len(hetero_expanded) > 0:
        # Merge with sequence data to get sequences
        hetero_expanded = pd.merge(
            hetero_expanded,
            sequence_df[['model_gene_id', 'sequence']],
            left_on='gene',
            right_on='model_gene_id',
            how='left'
        )
        
        hetero_expanded['enzyme_type'] = 'heteromeric'
        hetero_expanded['rate_value'] = hetero_expanded['SA_app_max']
        hetero_expanded['rate_units'] = 'μmol/mg/min'
        
        # Select relevant columns for heteromeric
        hetero_columns = [
            'sequence', 'gene', 'SMILES', 'rate_value', 'rate_units', 'enzyme_type',
            'rxn', 'condition_max', 'p_total_max', 'flux_value', 'subsystem',
            'eta_mean', 'eta_stdev', 'eta_min', 'eta_max', 'eta_cv',
            'SA_app_max', 'complex_genes', 'n_subunits', 'subunit_role',
            'subunit_molecular_weight', 'subunit_protein_ppm', 'complex_g_gdcw'
        ]
        hetero_columns = [col for col in hetero_columns if col in hetero_expanded.columns]
        hetero_view = hetero_expanded[hetero_columns].copy()
    else:
        hetero_view = pd.DataFrame()
    
    # Combine both views
    if len(hetero_view) > 0:
        unified_view = pd.concat([homo_view, hetero_view], ignore_index=True, sort=False)
    else:
        unified_view = homo_view.copy()
    
    # Sort by rate_value (descending)
    unified_view = unified_view.sort_values('rate_value', ascending=False).reset_index(drop=True)
    
    print(f"\nUnified view created:")
    print(f"  Homomeric sequences: {len(homo_view)}")
    print(f"  Heteromeric subunits: {len(hetero_view)}")
    print(f"  Total sequence-substrate pairs: {len(unified_view)}")
    
    # Check for sequences present in both homomeric and heteromeric
    if len(hetero_view) > 0:
        homo_genes = set(homo_view['gene'].dropna())
        hetero_genes = set(hetero_view['gene'].dropna())
        overlap = homo_genes & hetero_genes
        
        if overlap:
            print(f"\n  Note: {len(overlap)} genes appear in both homomeric and heteromeric forms")
            print(f"  Examples: {list(overlap)[:5]}")
    
    return unified_view


def convert_SA_to_kcat(unified_view: pd.DataFrame):
    """
    Add a converted kcat column for heteromeric enzymes to enable direct comparison
    with homomeric kcat values.
    
    Conversion: kcat [s⁻¹] = SA [μmol/mg/min] × MW [g/mol] / 1000 / 60
    
    For heteromeric subunits, uses the subunit molecular weight.
    
    Parameters:
        unified_view: pd.DataFrame
            Output from def create_unified_kcat_SA(kmax_results: pd.DataFrame, 

    
    Returns:
        pd.DataFrame: Same dataframe with added 'kcat_converted' column
    """
    
    print("Converting SA to kcat for heteromeric subunits...")
    
    df = unified_view.copy()
    
    # Initialize kcat_converted column
    df['kcat_converted'] = np.nan
    
    # For homomeric: kcat_converted = kcat_app_max (already in s⁻¹)
    if 'kcat_app_max' in df.columns:
        homo_mask = df['enzyme_type'] == 'homomeric'
        df.loc[homo_mask, 'kcat_converted'] = df.loc[homo_mask, 'kcat_app_max']
    
    # For heteromeric: convert SA to kcat using subunit MW
    hetero_mask = df['enzyme_type'] == 'heteromeric'
    
    if hetero_mask.any() and 'subunit_molecular_weight' in df.columns:
        # kcat [s⁻¹] = SA [μmol/mg/min] × MW [g/mol] / 1000 / 60
        # where MW in g/mol needs to be converted to mg/μmol
        # 1 g/mol = 1 mg/mmol = 1 mg/(1000 μmol) = 0.001 mg/μmol
        
        df.loc[hetero_mask, 'kcat_converted'] = (
            df.loc[hetero_mask, 'rate_value'] *  # SA in μmol/mg/min
            df.loc[hetero_mask, 'subunit_molecular_weight'] /  # MW in g/mol
            1000 /  # Convert to mg/μmol
            60      # Convert min to seconds
        )
        
        converted_count = df.loc[hetero_mask, 'kcat_converted'].notna().sum()
        print(f"  Converted SA to kcat for {converted_count} heteromeric subunits")
        
        if converted_count > 0:
            kcat_range = df.loc[hetero_mask & df['kcat_converted'].notna(), 'kcat_converted']
            print(f"  Converted kcat range: {kcat_range.min():.2e} to {kcat_range.max():.2e} s⁻¹")
    
    return df