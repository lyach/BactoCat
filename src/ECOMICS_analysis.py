"""
Module Description: ECOMICS_analysis.py

Purpose: Analyze ECOMICS fluxomics dataset.

This module provides functions to:
1. Select the most complete study from an ECOMICS dataset, based on non-zero fluxes, 
2. Modify the iML1515 GEM based on specified medium and stress conditions.
3. Create a dataframe that combines experimental fluxomics data with FBA and pFBA results, 
allowing for comparison and analysis of metabolic fluxes under different conditions.
"""

import pandas as pd
import cobra
from cobra.io import load_model
from cobra import flux_analysis


def select_most_complete_study(dataframe):
    """
    Select the most complete study from the ECOMICS dataset, based on the number of non-zero fluxes.
    
    Parameters:
        dataframe (pd.DataFrame): DataFrame with study conditions in first 4 columns and reaction fluxes in remaining columns
    
    Returns:
        pd.DataFrame: DataFrame with study conditions (Strain, MediumID, Medium, Stress), count of non-zero reactions, and all flux columns
    """
    # Extract the first 4 columns (study conditions)
    study_conditions = dataframe.iloc[:, :4].copy()
    
    # Extract reaction columns (from column 5 onwards)
    reaction_columns = dataframe.iloc[:, 4:]
    
    # Count non-zero values for each row (study)
    num_nonzero_rxns = (reaction_columns != 0.0).sum(axis=1)
    
    # Create result dataframe with study conditions, non-zero reaction count, and all flux columns
    result_df = study_conditions.copy()
    result_df['num_nonzero_rxns'] = num_nonzero_rxns
    result_df = pd.concat([result_df, reaction_columns], axis=1)
    
    return result_df


def apply_ecomics_condition(medium_id: str, stress: str):
    """
    Modify iML1515 GEM based on specified experimental conditions.
    
    Parameters:
        strain : str
            E. coli strain type. Options: 'W3110', 'MG1655', 'BW25113'
        medium_id : str
            Medium identifier. Options: 'MD066', 'MD120', 'MD004', 'MD121'
        medium : str
            Medium description. Options: 'synthetic+Glu', 'synthetic+Glu+KNO3', 'MOPS+Glu(0.4%)', 'M9+Glu'
        stress : str
            Stress condition. Options: 'none', 'NADH-limitation', 'ATP-limitation'
    
    Returns:
        modified_model: cobra.Model
            Modified GEM model
    """
    
    # Validate input parameters
    valid_medium_ids = {'MD066', 'MD120', 'MD004', 'MD121'}
    valid_stress = {'none', 'NADH-limitation', 'ATP-limitation'}
    
    if medium_id not in valid_medium_ids:
        raise ValueError(f"Invalid medium_id '{medium_id}'. Must be one of {valid_medium_ids}")
    if stress not in valid_stress:
        raise ValueError(f"Invalid stress '{stress}'. Must be one of {valid_stress}")
    
    # Load the iML1515 model
    model = cobra.io.load_model("iML1515")
    
    # Apply medium-specific modifications
    _apply_medium_conditions(model, medium_id)
    
    # Apply stress-specific modifications
    _apply_stress_conditions(model, stress)
    
    return model


def _apply_medium_conditions(model: cobra.Model, medium_id: str) -> None:
    """
    Apply medium-specific modifications to the model.
    
    Parameters:
        model : cobra.Model
            The genome-scale model to modify
        medium_id : str
            Medium identifier
    """
    
    if medium_id == 'MD066':
        # synthetic+Glu medium
        # TODO: add medium composition
        print("Applying synthetic+Glu medium")
        medium = {
            "EX_glc__D_e": -10
        }
        
    elif medium_id == 'MD004':
        # synthetic+Glu medium with higher glucose uptake
        # TODO: add medium composition
        print("Applying synthetic+Glu medium with higher glucose uptake")
        medium = {
            "EX_glc__D_e": -10
        }
    
    elif medium_id == 'MD120':
        # MOPS+Glu(0.4%) medium
        # TODO: add medium composition
        print("Applying MOPS+Glu(0.4%) medium")
        medium = {
            "EX_glc__D_e": -10
        }
    
    elif medium_id == 'MD121':
        # M9+Glu medium
        # NH4 and SO4 bounds were discarded, as they arent't growth limiting
        # Biomass flux and overall correlation improves by doing this
        print("Applying M9+Glu medium")
        medium = {
            "EX_glc__D_e": -10,
            #"EX_nh4_e": -5.229, 
            "EX_so4_e": -1.699,
            "EX_o2_e": -14.49,
            #"EX_co2_e": 16.22,
            #"EX_h2o_e": -6.96
        }
    
    # Apply the medium by setting lower bounds
    for rxn_id, uptake in medium.items():
        try:
            rxn = model.reactions.get_by_id(rxn_id)
            rxn.lower_bound = uptake
            print(f"Set {rxn_id} lower bound to {uptake}")
        except KeyError:
            print(f"Reaction {rxn_id} not found in model.")


def _apply_stress_conditions(model: cobra.Model, stress: str) -> None:
    """
    Apply stress-specific modifications to the model.
    
    Parameters:
        model : cobra.Model
            The genome-scale model to modify
        stress : str
            Stress condition
    """
    
    if stress == 'NADH-limitation':
        try:
            # Limit NADH dehydrogenase
            nadh_dehyd = model.reactions.get_by_id("NADH16pp")
            nadh_dehyd.upper_bound = 3
            print("Applied NADH dehydrogenase limitation")
            
            # Add NADH drain reaction
            nadh = model.metabolites.nadh_c
            dm_nadh = cobra.Reaction("DM_nadh_c")
            dm_nadh.name = "NADH drain reaction"
            dm_nadh.add_metabolites({nadh: -1})  # drains 1 NADH → nothing
            dm_nadh.lower_bound = dm_nadh.upper_bound = 0.1
            model.add_reactions([dm_nadh])
            print("Added NADH drain reaction")
            
        except KeyError as e:
            print(f"Error applying NADH limitation: {e}")
    
    elif stress == 'ATP-limitation':
        try:
            # Increase ATP maintenance requirement
            atp_maint = model.reactions.get_by_id("ATPM")
            atp_maint.lower_bound = 20  # mmol · 1/gDW · 1/h
            atp_maint.upper_bound = 20
            print("Increased ATP maintenance requirement")
            
            # Limit ATP synthase
            atp_synth = model.reactions.get_by_id("ATPS4rpp")
            atp_synth.upper_bound = 5
            print("Limited ATP synthase capacity")
            
        except KeyError as e:
            print(f"Error applying ATP limitation: {e}")
    
    elif stress == 'none':
        # No stress modifications
        pass
    
    
def create_fluxomics_dataframe(method: str, exp_fluxes: str, 
                               modified_gem: cobra.Model):
    """
    Create a dataframe with experimental, FBA and pFBA fluxomics results.
    
    Parameters:
        method: str
            Method for the flux simulations: 'FBA' or 'pFBA'
        exp_fluxes: str
            Path to the 'fluxomics_{condition_string}.csv' file
        modified_gem: cobra.Model
            Modified GEM from apply_ecomics_condition module
    
    Returns:
        pd.DataFrame: DataFrame with columns: 'rxn_id', 'exp_reaction', 'exp_flux', 'FBA_flux', 'pFBA_flux'
    """
    
    # Load experimental fluxomics data
    exp_fluxes_df = pd.read_csv(exp_fluxes)
    
    # Remove rows with missing experimental reaction names
    exp_fluxes_df = exp_fluxes_df.dropna(subset=['exp_reaction'])

    # Run FBA and pFBA with status checks
    if method == 'FBA':
        solution = modified_gem.optimize()
        print(f"FBA status: {solution.status}")
        if solution.status != 'optimal':
            print("FBA optimization did not succeed.")
            return None, None, None
    elif method == 'pFBA':
        solution = flux_analysis.pfba(modified_gem)
        print(f"pFBA status: {solution.status}")
        if solution.status != 'optimal':
            print("pFBA optimization did not succeed.")
            return None, None, None
    else:
        raise ValueError(f"Invalid method '{method}'. Must be 'FBA' or 'pFBA'.")
    
    # Extract fluxes and create DataFrames
    method_flux = f'{method}_flux'
    fluxomics_df = pd.DataFrame.from_dict(solution.fluxes.to_dict(), orient='index', columns=[method_flux])

    # Set index name
    fluxomics_df.index.name = 'rxn_id'
    
    # Process experimental data for matching
    # Extract base reaction names by removing prefixes like 'R0091_'
    exp_fluxes_df['base_reaction'] = exp_fluxes_df['exp_reaction'].str.replace(r'^R\d+_', '', regex=True)
    
    # Create a mapping from base reaction names to experimental data
    exp_mapping = {}
    for _, row in exp_fluxes_df.iterrows():
        base_rxn = row['base_reaction']
        # Skip if base_rxn is NaN or empty
        if pd.isna(base_rxn) or base_rxn == '':
            continue
            
        if base_rxn not in exp_mapping:
            exp_mapping[base_rxn] = []
        exp_mapping[base_rxn].append({
            'exp_reaction': row['exp_reaction'],
            'exp_flux': row['exp_flux']/10 # scale the experimental fluxes column by the glucose uptake rate (-10 mM/h)
        })
    
    # Track matched experimental reactions
    matched_exp_reactions = set()
    
    # Add experimental data to fluxomics_df
    exp_reactions = []
    exp_fluxes = []
    
    for rxn_id in fluxomics_df.index:
        if rxn_id in exp_mapping:
            # Find the experimental data with the closest match to the method flux, ignoring sign
            method_flux_value = abs(fluxomics_df.at[rxn_id, method_flux])
            closest_match = min(exp_mapping[rxn_id], key=lambda x: abs(abs(x['exp_flux']) - method_flux_value))
            
            exp_reactions.append(closest_match['exp_reaction'])
            exp_fluxes.append(closest_match['exp_flux'])
            
            # Track matched reactions
            matched_exp_reactions.add(rxn_id)
            
            # Print info about multiple matches
            if len(exp_mapping[rxn_id]) > 1:
                print(f"Multiple experimental reactions found for {rxn_id}: {[x['exp_reaction'] for x in exp_mapping[rxn_id]]}, used {closest_match['exp_reaction']}")
        else:
            exp_reactions.append(None)
            exp_fluxes.append(None)
    
    # Find unmatched experimental reactions (filter out NaN values)
    all_exp_base_reactions = set(exp_fluxes_df['base_reaction'].dropna().unique())
    unmatched_exp_reactions = all_exp_base_reactions - matched_exp_reactions
    
    # Print unmatched experimental reactions
    if unmatched_exp_reactions:
        print(f"\nUnmatched experimental reactions ({len(unmatched_exp_reactions)}):")
        # Filter out any remaining NaN values before sorting
        valid_unmatched = [rxn for rxn in unmatched_exp_reactions if pd.notna(rxn) and rxn != '']
        for base_rxn in sorted(valid_unmatched):
            # Get all experimental reaction names for this base reaction
            exp_rxn_names = exp_fluxes_df[exp_fluxes_df['base_reaction'] == base_rxn]['exp_reaction'].tolist()
            print(f"  {base_rxn}: {exp_rxn_names}")
    else:
        print("\nAll experimental reactions were matched!")
    
    print(f"\nMatching summary:")
    print(f"  Total experimental reactions: {len(exp_fluxes_df)}")
    print(f"  Unique base reaction names: {len(all_exp_base_reactions)}")
    print(f"  Matched reactions: {len(matched_exp_reactions)}")
    print(f"  Unmatched reactions: {len(unmatched_exp_reactions)}")
    
    # Add experimental columns to fluxomics_df
    fluxomics_df['exp_reaction'] = exp_reactions
    fluxomics_df['exp_flux'] = exp_fluxes
    
    # Reset index to make rxn_id a column
    fluxomics_df = fluxomics_df.reset_index()
    
    # Reorder columns
    fluxomics_df = fluxomics_df[['rxn_id', 'exp_reaction', 'exp_flux', method_flux]]
    
    return fluxomics_df


def correct_reversible_reactions(fluxomics_df: pd.DataFrame,
                                 method: str) -> pd.DataFrame:
    """
    Correct reversible reactions in the fluxomics dataframe.
    
    Parameters:
        fluxomics_df: pd.DataFrame
            DataFrame with experimental, FBA and pFBA fluxomics results
        method: str
            Method for the flux simulations: 'FBA_flux' or 'pFBA_flux'
    
    Returns:
        pd.DataFrame: DataFrame with corrected experimental flux sign
    """
    
    # Create a copy of the dataframe to avoid modifying the original
    corrected_df = fluxomics_df.copy()
    
    # Identify reversible reactions
    reversible_reactions = fluxomics_df[fluxomics_df[method] < 0]
    
    # Correct fluxes for reversible reactions
    for rxn_id, flux in reversible_reactions.iterrows():
        if flux[method] < 0 and corrected_df.at[rxn_id, 'exp_flux'] > 0:
            corrected_df.at[rxn_id, 'exp_flux'] = -corrected_df.at[rxn_id, 'exp_flux']
    
    return corrected_df
    