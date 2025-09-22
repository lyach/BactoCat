"""
Module Description: kapp_builder.py

Purpose: 
Script for building the kapp dataframe.

Overview: 
This module provides functions to:


"""

import cobra
import pandas as pd
from cobra import flux_analysis
from itertools import product
from Bio.SeqUtils.ProtParam import ProteinAnalysis

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
        
        # == Fluxomics info ==
        # Merge with specific flux condition
        flux_subset = fluxomics_df[['rxn_id', flux_col]].copy()
        flux_subset.rename(columns={flux_col: 'flux_value'}, inplace=True)
        
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
            'C(=O)O',      # formic acid
            'CO',          # methanol
            'CCO',         # ethanol
            'CC(=O)O',     # acetic acid
            'C',           # methane
            'CC',          # ethane
            'CCC',         # propane
            'N',           # ammonia (as N)
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

def evaluate_kapp_homomeric(kapp_results: dict, threshold: float = 1e6):
    """
    Evaluate kapp for homomeric enzymes by filtering out unrealistic high values.
    
    Parameters:
        kapp_results: dict
            Nested dictionary with structure {condition: {p_total: dataframe}}
        threshold: float
            Threshold for filtering kcat_app values (default: 1e6 s⁻¹)
            Rows with kcat_app > threshold will be removed
    
    Returns:
        kapp_filtered_results: dict
            Same structure as input but with filtered dataframes
    """
    
    print(f"Filtering kcat_app values above threshold: {threshold:.0e} s⁻¹")
    
    # Initialize output dictionary
    kapp_filtered_results = {}
    
    # Track filtering statistics
    total_original_rows = 0
    total_filtered_rows = 0
    total_removed_rows = 0
    
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
            
            # Filter out rows where kcat_app exceeds threshold
            # Keep rows where kcat_app is NaN, <= threshold, or missing
            mask = (df_filtered['kcat_app'].isna()) | (df_filtered['kcat_app'] <= threshold)
            df_filtered = df_filtered[mask]
            
            filtered_count = len(df_filtered)
            removed_count = original_count - filtered_count
            
            print(f"    Filtered rows: {filtered_count}")
            print(f"    Removed rows: {removed_count}")
            
            if removed_count > 0:
                # Show some statistics about removed values
                removed_values = df[~mask]['kcat_app'].dropna()
                if len(removed_values) > 0:
                    print(f"    Removed kcat_app range: {removed_values.min():.2e} to {removed_values.max():.2e} s⁻¹")
            
            # Update statistics
            total_original_rows += original_count
            total_filtered_rows += filtered_count
            total_removed_rows += removed_count
            
            # Store the filtered dataframe
            kapp_filtered_results[condition_name][p_total_value] = df_filtered
    
    # Print summary statistics
    print("\nFiltering Summary:")
    print(f"Total original rows: {total_original_rows}")
    print(f"Total filtered rows: {total_filtered_rows}")
    print(f"Total removed rows: {total_removed_rows}")
    if total_original_rows > 0:
        removal_percentage = (total_removed_rows / total_original_rows) * 100
        print(f"Percentage removed: {removal_percentage:.1f}%")
    
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
        'gene', 'rxn', 'flux_value', 'protein_mmol_gdcw', 'subsystem' # Additional useful columns
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