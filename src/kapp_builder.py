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
from cobra.io import read_sbml_model
from loguru import logger
from tqdm import tqdm
from src.FVA_analysis.utils import cobra_to_fva_problem


def process_medium_df(medium_df: pd.DataFrame):
    """
    Process a medium dataframe to extract condition dictionaries for model bound modification.
    
    Parameters
    ----------
    medium_df : pd.DataFrame
        DataFrame where each row is a condition to simulate.
    
    Returns
    -------
    list of tuple
        List of (condition_id, medium_dict) tuples, where medium_dict maps
        reaction IDs to flux values for that condition.
    """
    # Columns to exclude from reaction mapping
    exclude_cols = {'condition_id', 'avg_growth'}
    
    # Get reaction columns
    rxn_cols = [col for col in medium_df.columns if col not in exclude_cols]
    
    conditions = []
    for _, row in medium_df.iterrows():
        condition_id = str(row['condition_id'])
        medium_dict = {rxn_id: row[rxn_id] for rxn_id in rxn_cols}
        conditions.append((condition_id, medium_dict))
    
    return conditions


def create_fluxomics_dataframe(
    flux_method: str,
    GEM: cobra.Model,
    carbon_uptake: list = None,
    oxygen_uptake: list = None,
    carbon_exchange_rxn: str = "EX_glc__D_e",
    oxygen_exchange_rxn: str = "EX_o2_e",
    medium_df: pd.DataFrame = None,
    medium_upper_bound: bool = False,
):
    """
    Create a dataframe with FBA or pFBA fluxomics for all conditions.
    
    Supports two modes:
    1. Carbon/Oxygen mode: Uses combinations of carbon_uptake and oxygen_uptake rates
    2. Medium DataFrame mode: Uses conditions from medium_df (takes priority if provided)
    
    Parameters
    ----------
    flux_method : str
        Method for the flux simulations: 'FBA' or 'pFBA'
    GEM : cobra.Model
        The GEM model to perform flux analysis on
    carbon_uptake : list, optional
        List of carbon uptake rates to test (mmol/gDW/h)
    oxygen_uptake : list, optional
        List of oxygen uptake rates to test (mmol/gDW/h)
    carbon_exchange_rxn : str, optional
        Reaction ID for carbon exchange (default: 'EX_glc__D_e')
    oxygen_exchange_rxn : str, optional
        Reaction ID for oxygen exchange (default: 'EX_o2_e')
    medium_df : pd.DataFrame, optional
        DataFrame where each row is a condition to simulate
    medium_upper_bound : bool, optional
        If True, sets both lower and upper bounds when using medium_df mode
    
    Returns
    -------
    pd.DataFrame
        DataFrame with columns: 'rxn_id', 'flux_cond1', 'flux_cond2', ..., 'flux_condN'
    """
    # Initialize results dictionary
    flux_results = {}
    
    # Get all reaction IDs for the dataframe
    rxn_ids = [rxn.id for rxn in GEM.reactions]
    
    # Create conditions if medium_df is provided
    if medium_df is not None:
        logger.info("Using given medium mode for flux simulations")
        conditions = process_medium_df(medium_df)
        
        for condition_id, medium_dict in tqdm(conditions, desc="Flux conditions", unit="cond"):            
            # Create a copy of the model to avoid modifying the original
            model_copy = GEM.copy()
            
            # Apply medium conditions using modify_reaction_bounds
            modify_reaction_bounds(model_copy, medium_dict, medium_upper_bound=medium_upper_bound, verbose=True)
            
            # Run FBA or pFBA
            if flux_method == 'FBA':
                solution = model_copy.optimize()
            elif flux_method == 'pFBA':
                solution = flux_analysis.pfba(model_copy)
                objective_value = solution.fluxes['BIOMASS_Ec_iML1515_core_75p37M']
                logger.debug(f"Condition {condition_id} completed, growth rate = {objective_value:.4f} 1/h")
            else:
                raise ValueError(f"Invalid method '{flux_method}'. Must be 'FBA' or 'pFBA'.")
            
            # Store results
            col_name = f'flux_{condition_id}'
            if solution.status == 'optimal':
                flux_results[col_name] = [solution.fluxes[rxn_id] for rxn_id in rxn_ids]
            else:
                logger.warning(f"Condition {condition_id} optimization failed with status: {solution.status}")
                flux_results[col_name] = [0.0] * len(rxn_ids)
        
        num_conditions = len(conditions)
    
    # If carbon_uptake and oxygen_uptake are provided
    elif carbon_uptake is not None and oxygen_uptake is not None:
        logger.info("Using given carbon/oxygen combinations for flux simulations")
        uptake_combinations = list(product(carbon_uptake, oxygen_uptake))
        
        for i, (carbon_rate, oxygen_rate) in tqdm(enumerate(uptake_combinations, 1), total=len(uptake_combinations), desc="Flux conditions", unit="cond"):
            logger.debug(f"Processing condition {i}: Carbon={carbon_rate}, Oxygen={oxygen_rate}")
            
            # Create a copy of the model to avoid modifying the original
            model_copy = GEM.copy()
            
            # Set carbon uptake rate
            try:
                carbon_rxn = model_copy.reactions.get_by_id(carbon_exchange_rxn)
                carbon_rxn.lower_bound = -abs(carbon_rate)  # Negative for uptake
            except KeyError:
                logger.warning(f"Carbon uptake reaction '{carbon_exchange_rxn}' not found. Skipping carbon constraint.")
            
            # Set oxygen uptake rate
            try:
                oxygen_rxn = model_copy.reactions.get_by_id(oxygen_exchange_rxn)
                oxygen_rxn.lower_bound = -abs(oxygen_rate)  # Negative for uptake
            except KeyError:
                logger.warning(f"Oxygen uptake reaction '{oxygen_exchange_rxn}' not found. Skipping oxygen constraint.")
            
            # Run FBA or pFBA
            if flux_method == 'FBA':
                solution = model_copy.optimize()
                logger.debug(f"Condition {i} completed, growth rate = {solution.objective_value:.4f} 1/h")
            elif flux_method == 'pFBA':
                solution = flux_analysis.pfba(model_copy)
                objective_value = solution.fluxes['BIOMASS_Ec_iML1515_core_75p37M']
                logger.debug(f"Condition {i} completed, growth rate = {objective_value:.4f} 1/h")
            else:
                raise ValueError(f"Invalid method '{flux_method}'. Must be 'FBA' or 'pFBA'.")
            
            # Store results
            if solution.status == 'optimal':
                flux_results[f'flux_cond{i}'] = [solution.fluxes[rxn_id] for rxn_id in rxn_ids]
            else:
                logger.warning(f"Condition {i} optimization failed with status: {solution.status}")
                flux_results[f'flux_cond{i}'] = [0.0] * len(rxn_ids)
        
        num_conditions = len(uptake_combinations)
        
    else:
        raise ValueError(
            "Either 'medium_df', or both 'carbon_uptake' and 'oxygen_uptake' must be provided."
        )
    
    # Create the output dataframe
    fluxomics_df = pd.DataFrame({'rxn_id': rxn_ids})
    
    # Add flux columns for each condition
    for condition, fluxes in flux_results.items():
        fluxomics_df[condition] = fluxes
    
    logger.info(f"Fluxomics dataframe created with {num_conditions} conditions")
    return fluxomics_df

def modify_reaction_bounds(model, medium, medium_upper_bound=False, verbose=True):
    """
    Modify reaction bounds in a COBRA model based on medium conditions.
    
    Uses the COBRApy model.medium property, which:
      1. Closes ALL exchange reactions (lower_bound = 0)
      2. Opens only the ones specified with the given uptake rate
    
    Free metabolites (water, protons, CO2) are kept unconstrained
    regardless of what the medium dict specifies for them.
    
    Parameters
    ----------
    model : cobra.Model
        The COBRA model to modify in-place
    medium : dict
        Dictionary mapping exchange reaction IDs to uptake rates
        (positive values, as expected by model.medium)
    medium_upper_bound : bool, optional
        If True, also fixes the upper bound so the flux is locked
        at the specified uptake value
    verbose : bool, optional
        If True, prints information about modified reactions
    
    Returns
    -------
    None
        Modifies the model in-place
    """
    if medium is None:
        return
    
    FREE_METABOLITES = {
        'EX_h2o_e',
        'EX_h_e',
        'EX_co2_e',
    }
    DEFAULT_FREE_BOUND = 1000.0
    
    model_medium = model.medium
    
    for k in model_medium:
        model_medium[k] = 0
    
    for rxn_id, flux_value in medium.items():
        if rxn_id in FREE_METABOLITES:
            model_medium[rxn_id] = DEFAULT_FREE_BOUND
            if verbose:
                logger.debug(f"  {rxn_id}: free metabolite, left unconstrained")
            continue
        
        if rxn_id in model_medium:
            model_medium[rxn_id] = abs(float(flux_value))
        else:
            logger.warning(f"  Reaction '{rxn_id}' not in model exchanges, skipping")
    
    for rxn_id in FREE_METABOLITES:
        if rxn_id not in model_medium:
            try:
                model.reactions.get_by_id(rxn_id)
                model_medium[rxn_id] = DEFAULT_FREE_BOUND
            except KeyError:
                pass
    
    model.medium = model_medium
    
    if verbose:
        for rxn_id in medium:
            if rxn_id in FREE_METABOLITES:
                continue
            try:
                rxn = model.reactions.get_by_id(rxn_id)
                logger.debug(f"  Set {rxn_id}: lower={rxn.lower_bound}, upper={rxn.upper_bound}")
            except KeyError:
                pass
    
    if medium_upper_bound:
        for rxn_id, flux_value in medium.items():
            if rxn_id in FREE_METABOLITES:
                continue
            try:
                rxn = model.reactions.get_by_id(rxn_id)
                rxn.upper_bound = rxn.lower_bound
                if verbose:
                    logger.debug(f"  Fixed {rxn_id}: lower={rxn.lower_bound}, upper={rxn.upper_bound}")
            except KeyError:
                pass


def create_FVA_dataframe(
    GEM_path: str,
    carbon_uptake: list = None,
    oxygen_uptake: list = None,
    mu_fraction: float = 0.9,
    solver: str = 'cplex',
    carbon_exchange_rxn: str = "EX_glc__D_e",
    oxygen_exchange_rxn: str = "EX_o2_e",
    medium_df: pd.DataFrame = None,
    medium_upper_bound: bool = False,
):
    """
    Run FVA for all conditions. Supports two modes:
    1. Carbon/Oxygen mode: Uses combinations of carbon_uptake and oxygen_uptake rates
    2. Medium DataFrame mode: Uses conditions from medium_df (takes priority if provided)
    
    Parameters
    ----------
    GEM_path : str
        Path to the SBML model file (XML).
    carbon_uptake : list, optional
        List of carbon uptake rates to test (mmol/gDW/h)
    oxygen_uptake : list, optional
        List of oxygen uptake rates to test (mmol/gDW/h)
    mu_fraction : float, optional
        Fraction of optimal growth rate for FVA (default = 0.9).
    solver : str, optional
        Solver to use ('cplex' or 'gurobi'). Default is 'cplex'.
    carbon_exchange_rxn : str, optional
        Reaction ID for carbon exchange (default: 'EX_glc__D_e').
    oxygen_exchange_rxn : str, optional
        Reaction ID for oxygen exchange (default: 'EX_o2_e').
    medium_df : pd.DataFrame, optional
        DataFrame where each row is a condition to simulate
    medium_upper_bound : bool, optional
        If True, sets both lower and upper bounds when using medium_df mode
    
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

    # Load base model once
    base_model = read_sbml_model(GEM_path)
    rxn_ids = [rxn.id for rxn in base_model.reactions]
    
    # Initialize dictionaries for lower/upper bounds
    FVA_lower_results = {}
    FVA_upper_results = {}
    
    # If medium_df is provided
    if medium_df is not None:
        logger.info("Using given medium mode for FVA simulations")
        conditions = process_medium_df(medium_df)
        
        for condition_id, medium_dict in tqdm(conditions, desc="FVA conditions", unit="cond"):            
            # Copy model
            model_copy = base_model.copy()
            
            # Apply medium conditions using modify_reaction_bounds
            modify_reaction_bounds(model_copy, medium_dict, medium_upper_bound=medium_upper_bound, verbose=True)
            
            # Optimize
            # TO DO - is it necessary to do prev optimization for FVA?
            solution = model_copy.optimize()
            print(f"Condition {condition_id}: Growth Rate = {solution.objective_value:.4f} 1/h")
            if solution.status != 'optimal':
                logger.warning(f"Optimization failed at condition {condition_id} with status: {solution.status}, filling with NaNs.")
                FVA_lower_results[f'FVA_lower_{condition_id}'] = [float('nan')] * len(rxn_ids)
                FVA_upper_results[f'FVA_upper_{condition_id}'] = [float('nan')] * len(rxn_ids)
                continue
            
            # Build FVA problem
            problem = cobra_to_fva_problem(model_copy, mu=mu_fraction)
            
            # Run FVA
            fva_results = fva_solve_faster(problem)
            fva_df = pd.DataFrame({
                'rxn_id': [rxn.id for rxn in model_copy.reactions],
                'FVA_lower': fva_results.lower_bound,
                'FVA_upper': fva_results.upper_bound
            }) 
            fva_df = fva_df.set_index('rxn_id').reindex(rxn_ids).reset_index()
            
            # Store FVA lower/upper bounds
            FVA_lower_results[f'FVA_lower_{condition_id}'] = fva_df['FVA_lower'].values
            FVA_upper_results[f'FVA_upper_{condition_id}'] = fva_df['FVA_upper'].values
            
        
        num_conditions = len(conditions)
        
    # If carbon_uptake and oxygen_uptake are provided
    elif carbon_uptake is not None and oxygen_uptake is not None:
        logger.info("Using given carbon/oxygen combinations for FVA simulations")
        uptake_combinations = list(product(carbon_uptake, oxygen_uptake))
        
        for i, (carbon_rate, oxygen_rate) in tqdm(enumerate(uptake_combinations, 1), total=len(uptake_combinations), desc="FVA conditions", unit="cond"):            
            # Copy model
            model_copy = base_model.copy()
            
            # Set medium
            try:
                carbon_rxn = model_copy.reactions.get_by_id(carbon_exchange_rxn)
                carbon_rxn.lower_bound = -abs(carbon_rate)
            except KeyError:
                logger.warning(f"Carbon uptake reaction '{carbon_exchange_rxn}' not found. Skipping.")
            try:
                oxygen_rxn = model_copy.reactions.get_by_id(oxygen_exchange_rxn)
                oxygen_rxn.lower_bound = -abs(oxygen_rate)
            except KeyError:
                logger.warning(f"Oxygen uptake reaction '{oxygen_exchange_rxn}' not found. Skipping.")
            
            # Optimize and get optimal mu
            solution = model_copy.optimize()
            if solution.status != 'optimal':
                logger.warning(f"Optimization failed at condition {i} with status: {solution.status}, filling with NaNs.")
                FVA_lower_results[f'FVA_lower_cond{i}'] = [float('nan')] * len(rxn_ids)
                FVA_upper_results[f'FVA_upper_cond{i}'] = [float('nan')] * len(rxn_ids)
                continue
            
            # Build FVA problem
            problem = cobra_to_fva_problem(model_copy, mu=mu_fraction)
            
            # Run FVA
            fva_results = fva_solve_faster(problem)
            fva_df = pd.DataFrame({
                'rxn_id': [rxn.id for rxn in model_copy.reactions],
                'FVA_lower': fva_results.lower_bound,
                'FVA_upper': fva_results.upper_bound
            }) 
            fva_df = fva_df.set_index('rxn_id').reindex(rxn_ids).reset_index()
            
            # Store FVA lower/upper bounds
            FVA_lower_results[f'FVA_lower_cond{i}'] = fva_df['FVA_lower'].values
            FVA_upper_results[f'FVA_upper_cond{i}'] = fva_df['FVA_upper'].values
                    
        num_conditions = len(uptake_combinations)
        
    else:
        raise ValueError(
            "Either 'medium_df' must be provided, or both 'carbon_uptake' and 'oxygen_uptake' must be provided."
        )
    
    # Build the output dataframe
    all_data = {'rxn_id': rxn_ids, **FVA_lower_results, **FVA_upper_results}
    fva_combined = pd.DataFrame(all_data)
    
    logger.info(f"FVA dataframe created with {num_conditions} conditions")
    return fva_combined


def FVA_integration(fluxomics_df: pd.DataFrame, fva_df: pd.DataFrame, filter: bool = False):
    """
    Check fluxes against FVA bounds for all conditions and optionally filter out reactions with violations.

    Returns:
        filtered_fluxomics_df, violations_df
    """
    merged_df = fluxomics_df.merge(fva_df, on='rxn_id', how='left')
    
    # Identify flux columns
    flux_cols = [col for col in merged_df.columns if col.startswith('flux_')]
    
    violations = []
    
    for col in flux_cols:
        # Map flux column to corresponding FVA lower/upper columns
        if col.startswith('flux_cond'):
            cond_suffix = col.replace('flux_', '')
        else:
            cond_suffix = col.replace('flux_', '')  # if using medium_df
        
        lower_col = f'FVA_lower_{cond_suffix}'
        upper_col = f'FVA_upper_{cond_suffix}'
        
        # Skip if FVA columns don't exist for this flux column
        if lower_col not in merged_df.columns or upper_col not in merged_df.columns:
            continue
        
        # Check lower bound violations - with small tolerance
        below_mask = merged_df[col] < (merged_df[lower_col]- 1e-6)
        above_mask = merged_df[col] > (merged_df[upper_col]+ 1e-6)
        
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
    
    # Violation report
    if not violations_df.empty:
        unique_rxns = violations_df['rxn_id'].nunique()
        total_violations = len(violations_df)
        logger.debug(f"Total violation instances (Cell count): {total_violations}")
        logger.debug(f"Unique reactions affected (Row count): {unique_rxns}")
        logger.debug(f"Violations by Condition: {(violations_df.groupby('condition').size().to_string())}")
    else:
        logger.debug("No FVA violations detected.")
    
    
    if filter and not violations_df.empty:
        violating_rxns = violations_df['rxn_id'].unique()
        before = len(merged_df)
        merged_df = merged_df[~merged_df['rxn_id'].isin(violating_rxns)].copy()
        after = len(merged_df)
        logger.debug(f"Filtered out {before - after} reactions with violations.")
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
    flux_columns = [col for col in fluxomics_df.columns if col.startswith('flux_')]
        
    # Process each flux condition
    for flux_col in flux_columns:
        logger.debug(f"Processing {flux_col}...")
        logger.debug(f" Rows before filtering: {len(enzymes_df)}")
        
        # Create a copy of enzymes_df for this condition
        condition_df = enzymes_df.copy()
        
        # Extract condition suffix
        cond_suffix = flux_col.replace("flux_", "")
        lower_col = f"FVA_lower_{cond_suffix}"
        upper_col = f"FVA_upper_{cond_suffix}"


        # == Fluxomics info ==
        # Merge with specific flux condition
        if lower_col not in fluxomics_df.columns or upper_col not in fluxomics_df.columns:
            logger.warning(f"Warning: Missing FVA columns for {flux_col} — skipping FVA merge.")
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
        
        logger.debug(f"  {flux_col}: {len(condition_df)} rows after filtering")
    
    logger.info(f"Created enzyme info dataframes for {len(enzyme_info_dfs)} conditions")
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


def process_enzyme_protein_mapping(enzyme_info_dfs: dict, paxdb_path: str, p_total: float):
    """
    Apply PaxDB protein mapping across all enzyme info dataframes.
    
    Parameters:
        enzyme_info_dfs: dict
            Dictionary with condition names as keys and enzyme dataframes as values
        paxdb_path: str
            Path to PaxDB TSV file
        p_total: float
            Total protein content value (g protein / g DCW)
    
    Returns:
        enzyme_protein_info_dfs: dict
            Dictionary structure: {condition_name: mapped_dataframe}
    """
    
    # Load PaxDB data
    logger.info(f"Loading PaxDB data from: {paxdb_path}")
    paxdb_df = pd.read_csv(paxdb_path, sep="\t", comment="#", header=None, names=["gene_name", "string_external_id", "abundance"])
    logger.debug(f"PaxDB data loaded: {len(paxdb_df)} rows")
    
    # Initialize output dictionary
    enzyme_protein_info_dfs = {}
    
    logger.info(f"Processing {len(enzyme_info_dfs)} conditions with p_total={p_total}")
    
    for condition_name, enzyme_df in enzyme_info_dfs.items():
        # Apply map_paxdb_to_gene function
        try:
            mapped_df = map_paxdb_to_gene(
                paxdb_df=paxdb_df,
                df_enzymes=enzyme_df,
                p_total=p_total
            )
            
            # Store the result
            enzyme_protein_info_dfs[condition_name] = mapped_df
            logger.debug(f"  {condition_name}: {len(mapped_df)} rows, {mapped_df['protein_ppm'].notna().sum()} with protein data")
            
        except Exception as e:
            logger.error(f"Error processing {condition_name}: {str(e)}")
            # Store None for failed conditions
            enzyme_protein_info_dfs[condition_name] = None
    
    return enzyme_protein_info_dfs



def calculate_kapp_homomeric(enzyme_protein_info_dfs: dict):
    """
    Calculate kapp for homomeric enzymes for each condition.
    
    Parameters:
        enzyme_protein_info_dfs: dict
            Dictionary with structure {condition: dataframe}
    
    Returns:
        dict: Same structure as input but with added 'kcat_app' column in each dataframe
    """
    
    # Initialize output dictionary
    kapp_results = {}
    
    for condition_name, enzyme_df in enzyme_protein_info_dfs.items():
        # Skip if dataframe is None (failed processing)
        if enzyme_df is None:
            logger.debug(f"  {condition_name}: Skipping - no data available")
            kapp_results[condition_name] = None
            continue
        
        # Work with a copy to avoid modifying original
        df_copy = enzyme_df.copy()
        initial_rows = len(df_copy)
        
        # Keep only homomeric enzymes
        df_copy = df_copy[df_copy['gpr_class'] == 'simple']
        
        # Drop duplicate enzymes
        df_copy = df_copy.drop_duplicates(subset=["gene", "SMILES"])

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
        logger.debug(f"  {condition_name}: {initial_rows} → {len(df_copy)} rows (homomeric+dedup), {valid_kcat} valid kcat_app")
        
        # Store the processed dataframe
        kapp_results[condition_name] = df_copy
    
    logger.info("Completed kcat_app calculation for all conditions")
    return kapp_results

def evaluate_kapp_homomeric(kapp_results: dict, upper_threshold: float = 1e6, lower_threshold: float = 1e-5):
    """
    Evaluate kapp for homomeric enzymes by filtering out unrealistic high and low values.
    
    Parameters:
        kapp_results: dict
            Dictionary with structure {condition: dataframe}
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
    logger.info(f"Filtering kcat_app values outside range: {lower_threshold:.0e} to {upper_threshold:.0e} s⁻¹")
    
    # Initialize output dictionary
    kapp_filtered_results = {}
    
    # Track filtering statistics
    total_original_rows = 0
    total_filtered_rows = 0
    total_removed_high = 0
    total_removed_low = 0
    
    for condition_name, df in kapp_results.items():
        # Skip if dataframe is None
        if df is None:
            kapp_filtered_results[condition_name] = None
            continue
        
        # Work with a copy to avoid modifying original
        df_filtered = df.copy()
        original_count = len(df_filtered)
        
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
        
        # Update statistics
        total_original_rows += original_count
        total_filtered_rows += filtered_count
        total_removed_high += removed_high_count
        total_removed_low += removed_low_count
        
        # Store the filtered dataframe
        kapp_filtered_results[condition_name] = df_filtered
    
    # Log summary statistics
    total_removed_rows = total_removed_high + total_removed_low
    if total_original_rows > 0:
        removal_percentage = (total_removed_rows / total_original_rows) * 100
        logger.info(
            f"Filtering complete: {total_original_rows} → {total_filtered_rows} rows "
            f"(removed {total_removed_rows}: {total_removed_high} high, {total_removed_low} low, {removal_percentage:.1f}%)"
        )
    
    return kapp_filtered_results

def get_kmax_homomeric(kapp_results: dict):
    """
    Get the maximum kapp of each enzyme-substrate pair across all conditions.
    
    Parameters:
        kapp_results: dict
            Dictionary with structure {condition: dataframe}
    Returns:
        kmax_results: pd.DataFrame
            DataFrame with columns: ['sequence', 'SMILES', 'kcat_app_max', 'condition_max']
            containing the maximum kcat_app value for each enzyme-substrate pair
    """
    logger.info("Starting kmax analysis across all conditions...")
    
    # List to collect all dataframes with metadata
    all_dataframes = []
    
    for condition_name, df in kapp_results.items():
        # Skip None dataframes
        if df is None or len(df) == 0:
            continue
            
        # Add metadata columns to track source
        df_with_metadata = df.copy()
        df_with_metadata['source_condition'] = condition_name
        
        # Only keep rows with valid kcat_app values
        df_with_metadata = df_with_metadata[df_with_metadata['kcat_app'].notna()]
        
        if len(df_with_metadata) > 0:
            all_dataframes.append(df_with_metadata)
    
    if not all_dataframes:
        logger.warning("No valid data found across all conditions")
        return pd.DataFrame(columns=['sequence', 'SMILES', 'kcat_app_max', 'condition_max'])
    
    # Concatenate all dataframes
    combined_df = pd.concat(all_dataframes, ignore_index=True)
    logger.debug(f"Combined dataframe has {len(combined_df)} total entries")
    
    # Group by enzyme-substrate pair (sequence + SMILES) and find maximum kcat_app
    kmax_results = (
        combined_df.loc[combined_df.groupby(['sequence', 'SMILES'])['kcat_app'].idxmax()]
        .reset_index(drop=True)
    )
    
    # Select and rename relevant columns for output
    output_columns = [
        'sequence', 'SMILES', 'kcat_app', 'source_condition',
        'gene', 'rxn', 'flux_value', 'FVA_upper', 'FVA_lower', 'protein_mmol_gdcw', 'subsystem' 
    ]
    
    # Keep only columns that exist in the dataframe
    available_columns = [col for col in output_columns if col in kmax_results.columns]
    kmax_results = kmax_results[available_columns].copy()
    
    # Rename columns for clarity
    column_renames = {
        'kcat_app': 'kcat_app_max',
        'source_condition': 'condition_max'
    }
    kmax_results = kmax_results.rename(columns=column_renames)
    
    # Sort by kcat_app_max in descending order
    kmax_results = kmax_results.sort_values('kcat_app_max', ascending=False).reset_index(drop=True)
    
    logger.info(
        f"Found kmax for {len(kmax_results)} enzyme-substrate pairs "
        f"(range: {kmax_results['kcat_app_max'].min():.2e} to {kmax_results['kcat_app_max'].max():.2e} s⁻¹)"
    )
    
    return kmax_results


def get_eta(kapp_results: dict, kmax_results: pd.DataFrame):
    """
    Calculate eta (kapp/kmax) for each enzyme-substrate pair across all conditions.
    
    Parameters:
        kapp_results: dict
            Dictionary with structure {condition: dataframe}
        kmax_results: pd.DataFrame
            DataFrame with maximum kcat_app values for each enzyme-substrate pair
    
    Returns:
        tuple: (kapp_results_with_eta, kmax_results_with_variance)
            - kapp_results_with_eta: dict with same structure as input but with 'eta' column added
            - kmax_results_with_variance: DataFrame with added variance columns (eta_mean, eta_stdev, eta_min, eta_max, eta_cv)
    """
    logger.info("Calculating eta (kapp/kmax) for all conditions...")
    
    # Initialize output dictionary
    kapp_results_with_eta = {}
    
    # List to collect all eta values for variance calculation
    all_eta_values = []
    
    for condition_name, df in kapp_results.items():
        # Skip if dataframe is None
        if df is None:
            kapp_results_with_eta[condition_name] = None
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
        
        # Store the dataframe with eta
        kapp_results_with_eta[condition_name] = df_with_eta
        
        # Collect eta values for variance calculation
        eta_data = df_with_eta[['sequence', 'SMILES', 'eta']].copy()
        eta_data['source_condition'] = condition_name
        eta_data = eta_data[eta_data['eta'].notna()]  # Keep only valid eta values
        
        if len(eta_data) > 0:
            all_eta_values.append(eta_data)
    
    # Calculate variance metrics for each enzyme-substrate pair
    if not all_eta_values:
        logger.warning("No valid eta values found")
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
    
    # Log summary statistics
    logger.info(
        f"Eta calculation complete: {len(kmax_with_variance)} pairs, "
        f"mean η={variance_metrics['eta_mean'].mean():.3f}, CV={variance_metrics['eta_cv'].mean():.3f}"
    )
    
    return kapp_results_with_eta, kmax_with_variance
