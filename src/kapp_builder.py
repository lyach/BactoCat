"""
kapp_builder.py

Module with functions to build the kmax and kapp DataFrames.
"""

import cobra
import pandas as pd
from cobra import flux_analysis
from cobra.io import read_sbml_model
from loguru import logger
from tqdm import tqdm

from src.FVA_analysis.utils import cobra_to_fva_problem
from src.proteomics_mapper import calculate_molecular_weight
from src.utils import load_dataframe_if_path

# =============================================================================
# Constants
# =============================================================================

DEFAULT_FREE_BOUND = 1000.0

COFACTORS = {
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

GROWTH_REACTION = 'BIOMASS_Ec_iML1515_core_75p37M'

# =============================================================================
# Functions
# =============================================================================

def process_medium_df(medium_df: pd.DataFrame, 
                      include_growth: bool = True,
                      growth_key: str = 'avg_growth',
                      growth_reaction: str = GROWTH_REACTION):
    """
    Process a medium dataframe to extract condition dictionaries for model bound modification.
    
    Parameters
    ----------
    medium_df : pd.DataFrame
        DataFrame where each row is a condition to simulate.
    include_growth : bool, optional
        If True, include growth value as a key in the condition dictionary.
    growth_key : str
        The key in the medium dataframe that contains the growth value.
    growth_reaction : str
        The reaction ID of the growth reaction.
    Returns
    -------
    list of tuple
        List of (condition_id, medium_dict) tuples, where medium_dict maps
        reaction IDs to flux values for that condition.
    """
    exclude_cols = {'condition_id', growth_key}
    
    rxn_cols = [col for col in medium_df.columns if col not in exclude_cols]
    
    conditions = []
    for _, row in medium_df.iterrows():
        condition_id = str(row['condition_id'])
        medium_dict = {rxn_id: row[rxn_id] for rxn_id in rxn_cols}
        if include_growth and growth_reaction is not None and growth_key in medium_df.columns:
            medium_dict[growth_reaction] = row[growth_key]
        conditions.append((condition_id, medium_dict))
    
    return conditions


def create_fluxomics_dataframe(
    flux_method: str,
    GEM: cobra.Model,
    medium_df: pd.DataFrame,
    medium_upper_bound: bool = False,
    growth_reaction: str = GROWTH_REACTION,
    include_growth: bool = False,
    free_metabolites=None,
):
    """
    Create a dataframe with FBA or pFBA fluxomics for all conditions.

    Parameters
    ----------
    flux_method : str
        Method for the flux simulations: 'FBA' or 'pFBA'
    GEM : cobra.Model
        The GEM model to perform flux analysis on
    medium_df : pd.DataFrame
        DataFrame where each row is a condition to simulate. Must contain a
        ``condition_id`` column.
    medium_upper_bound : bool, optional
        If True, sets both lower and upper bounds when applying medium conditions.
    free_metabolites : list, optional
        Exchange reaction IDs to leave unconstrained.

    Returns
    -------
    pd.DataFrame
        DataFrame with columns: 'rxn_id', '<condition_id_1>', '<condition_id_2>', ...
    """
    if medium_df is None:
        raise ValueError("'medium_df' must be provided.")

    # Initialize results dictionary
    flux_results = {}

    # Get all reaction IDs for the dataframe
    rxn_ids = [rxn.id for rxn in GEM.reactions]

    logger.info("Using given medium mode for flux simulations")
    conditions = process_medium_df(medium_df, include_growth=include_growth, growth_reaction=growth_reaction)

    failed_conditions = []

    for condition_id, medium_dict in tqdm(conditions, desc="Flux conditions", unit="cond"):
        model_copy = GEM.copy()

        modify_reaction_bounds(model_copy, medium_dict, medium_upper_bound=medium_upper_bound,
                               growth_reaction=growth_reaction, verbose=True,
                               free_metabolites=free_metabolites)

        try:
            if flux_method == 'FBA':
                solution = model_copy.optimize()
            elif flux_method == 'pFBA':
                solution = flux_analysis.pfba(model_copy)
            else:
                raise ValueError(f"Invalid method '{flux_method}'. Must be 'FBA' or 'pFBA'.")
        except (cobra.exceptions.Infeasible, cobra.exceptions.OptimizationError) as e:
            logger.warning(f"Condition '{condition_id}' failed optimization: {e} — skipping.")
            failed_conditions.append(condition_id)
            continue

        if solution.status == 'optimal':
            flux_results[condition_id] = [solution.fluxes[rxn_id] for rxn_id in rxn_ids]
        else:
            logger.warning(f"Condition '{condition_id}' returned status '{solution.status}' — skipping.")
            failed_conditions.append(condition_id)

    if failed_conditions:
        logger.warning(
            f"{len(failed_conditions)}/{len(conditions)} conditions failed and were discarded: "
            f"{failed_conditions}"
        )

    if not flux_results:
        raise RuntimeError("All conditions failed optimization — cannot proceed.")

    fluxomics_df = pd.concat(
        [pd.DataFrame({'rxn_id': rxn_ids}), pd.DataFrame(flux_results, index=range(len(rxn_ids)))],
        axis=1,
    )

    logger.info(f"Fluxomics dataframe created with {len(flux_results)} successful conditions "
                f"(out of {len(conditions)} total)")
    return fluxomics_df


def modify_reaction_bounds(model, medium, medium_upper_bound=False,
                           growth_reaction=GROWTH_REACTION, verbose=True,
                           free_metabolites=None):
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
        (positive values, as expected by model.medium).
    medium_upper_bound : bool, optional
        If True, also fixes the upper bound so the flux is locked
        at the specified uptake value
    growth_reaction : str, optional
        Reaction ID of the biomass/growth reaction. If present in medium,
        its lower and upper bounds are fixed to the measured growth value.
    verbose : bool, optional
        If True, prints information about modified reactions
    free_metabolites : list or set, optional
        Exchange reaction IDs to leave unconstrained.
    Returns
    -------
    None
        Modifies the model in-place
    """
    if medium is None:
        return

    _free = set(free_metabolites) if free_metabolites else set()

    exchange_ids = {rxn.id for rxn in model.exchanges}
    
    model_medium = model.medium
    for k in model_medium:
        model_medium[k] = 0
    
    for rxn_id, flux_value in medium.items():
        if rxn_id == growth_reaction:
            continue
        
        if rxn_id in _free:
            model_medium[rxn_id] = DEFAULT_FREE_BOUND
            if verbose:
                logger.debug(f"  {rxn_id}: free metabolite, left unconstrained")
            continue
        
        if pd.isna(flux_value):
            continue
        if rxn_id in exchange_ids:
            model_medium[rxn_id] = abs(float(flux_value))
        else:
            logger.warning(f"  Reaction '{rxn_id}' not found in model, skipping")
    
    for rxn_id in _free:
        if rxn_id in exchange_ids:
            model_medium[rxn_id] = DEFAULT_FREE_BOUND
    
    model.medium = model_medium
    
    # Fix growth to measured value
    if growth_reaction is not None and growth_reaction in medium:
        try:
            rxn = model.reactions.get_by_id(growth_reaction)
            growth_value = float(medium[growth_reaction])
            rxn.lower_bound = growth_value
            rxn.upper_bound = growth_value
            if verbose:
                logger.debug(f"  Fixed growth {growth_reaction}: lower={rxn.lower_bound}, upper={rxn.upper_bound}")
        except KeyError:
            logger.warning(f"  Growth reaction '{growth_reaction}' not found in model, skipping")
    
    if verbose:
        for rxn_id in medium:
            if rxn_id in _free:
                continue
            try:
                rxn = model.reactions.get_by_id(rxn_id)
                logger.debug(f"  Set {rxn_id}: lower={rxn.lower_bound}, upper={rxn.upper_bound}")
            except KeyError:
                pass
    
    if medium_upper_bound:
        for rxn_id, flux_value in medium.items():
            if rxn_id in _free:
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
    medium_df: pd.DataFrame,
    mu_fraction: float = 0.9,
    solver: str = 'cplex',
    medium_upper_bound: bool = False,
    growth_reaction: str = GROWTH_REACTION,
    include_growth: bool = False,
    free_metabolites=None,
):
    """
    Run FVA for all conditions defined in ``medium_df``.

    Parameters
    ----------
    GEM_path : str
        Path to the SBML model file (XML).
    medium_df : pd.DataFrame
        DataFrame where each row is a condition to simulate. Must contain a
        ``condition_id`` column.
    mu_fraction : float, optional
        Fraction of optimal growth rate for FVA (default = 0.9).
    solver : str, optional
        Solver to use ('cplex' or 'gurobi'). Default is 'cplex'.
    medium_upper_bound : bool, optional
        If True, sets both lower and upper bounds when applying medium conditions.
    free_metabolites : list, optional
        Exchange reaction IDs to leave unconstrained.

    Returns
    -------
    pd.DataFrame
        Combined FVA dataframe with columns:
        ['rxn_id', 'FVA_lower_<cond>', 'FVA_upper_<cond>', ...]
    """
    if medium_df is None:
        raise ValueError("'medium_df' must be provided.")
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
    
    logger.info("Using given medium mode for FVA simulations")
    conditions = process_medium_df(medium_df, include_growth=include_growth, growth_reaction=growth_reaction)

    failed_conditions = []

    for condition_id, medium_dict in tqdm(conditions, desc="FVA conditions", unit="cond"):
        model_copy = base_model.copy()

        modify_reaction_bounds(model_copy, medium_dict, medium_upper_bound=medium_upper_bound,
                               growth_reaction=growth_reaction, verbose=True,
                               free_metabolites=free_metabolites)

        try:
            problem = cobra_to_fva_problem(model_copy, mu=mu_fraction)
            fva_results = fva_solve_faster(problem)
        except (cobra.exceptions.Infeasible, cobra.exceptions.OptimizationError, Exception) as e:
            logger.warning(f"FVA condition '{condition_id}' failed: {e} — skipping.")
            failed_conditions.append(condition_id)
            continue

        fva_df = pd.DataFrame({
            'rxn_id': [rxn.id for rxn in model_copy.reactions],
            'FVA_lower': fva_results.lower_bound,
            'FVA_upper': fva_results.upper_bound
        })
        fva_df = fva_df.set_index('rxn_id').reindex(rxn_ids).reset_index()

        FVA_lower_results[f'FVA_lower_{condition_id}'] = fva_df['FVA_lower'].values
        FVA_upper_results[f'FVA_upper_{condition_id}'] = fva_df['FVA_upper'].values

    if failed_conditions:
        logger.warning(
            f"FVA: {len(failed_conditions)}/{len(conditions)} conditions failed and were discarded: "
            f"{failed_conditions}"
        )

    if not FVA_lower_results:
        raise RuntimeError("All FVA conditions failed — cannot proceed.")

    all_data = {'rxn_id': rxn_ids, **FVA_lower_results, **FVA_upper_results}
    fva_combined = pd.DataFrame(all_data)

    logger.info(f"FVA dataframe created with {len(FVA_lower_results)} successful conditions "
                f"(out of {len(conditions)} total)")
    return fva_combined


def FVA_integration(fluxomics_df: pd.DataFrame, fva_df: pd.DataFrame, filter: bool = False):
    """
    Check fluxes against FVA bounds for all conditions and optionally filter out reactions with violations.

    Returns:
        filtered_fluxomics_df, violations_df
    """
    merged_df = fluxomics_df.merge(fva_df, on='rxn_id', how='left')
    
    # Identify condition columns
    flux_cols = [col for col in merged_df.columns if col != 'rxn_id' and not col.startswith('FVA_')]

    violations = []

    for col in flux_cols:
        lower_col = f'FVA_lower_{col}'
        upper_col = f'FVA_upper_{col}'
        
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


def create_enzyme_info_dataframe(enzymes_df, fluxomics_df, substrates_df, sequence_df, run_fva: bool = True):
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
        run_fva: bool
            If True, run flux variability analysis (default: True)
    
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
    
    # Get all condition columns — everything except rxn_id and FVA bound columns
    flux_columns = [col for col in fluxomics_df.columns if col != 'rxn_id' and not col.startswith('FVA_')]
        
    # Process each flux condition
    for flux_col in flux_columns:
        logger.debug(f"Processing {flux_col}...")
        logger.debug(f" Rows before filtering: {len(enzymes_df)}")
        
        # Create a copy of enzymes_df for this condition
        condition_df = enzymes_df.copy()
        
        lower_col = f"FVA_lower_{flux_col}"
        upper_col = f"FVA_upper_{flux_col}"


        # Fluxomics info
        # Merge with specific flux condition
        if lower_col not in fluxomics_df.columns or upper_col not in fluxomics_df.columns:
            if run_fva:
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
        
        # Get substrate info
        substrates_clean = substrates_df[['Reaction', 'SMILES', 'Direction']].copy()
        
        # Merge substrate info
        condition_df = pd.merge(condition_df, substrates_clean, left_on="rxn", right_on="Reaction", how="left")
        
        # Sequence info
        # Merge sequence info
        condition_df = pd.merge(condition_df, sequence_df, left_on="gene", right_on="model_gene_id", how="left")

        # Molecular weight (g/mol) from sequence
        #  required for the mass-fraction calculation of isoenzymes and complexes
        if 'sequence' in condition_df.columns:
            condition_df['molecular_weight'] = condition_df['sequence'].apply(calculate_molecular_weight)
        
        # Data cleaning and filtering
        # Drop rows with wrong direction-flux
        condition_df = condition_df[
            ((condition_df['Direction'] == 'forward') & (condition_df['flux_value'] >= 0)) |
            ((condition_df['Direction'] == 'reverse') & (condition_df['flux_value'] <= 0))
        ]
        
        # Drop rows with balancing species / cofactors (unlikely to be substrates)
        condition_df = condition_df[~condition_df['SMILES'].isin(COFACTORS)]
        
        # Remove rows with 0 flux 
        condition_df = condition_df[condition_df['flux_value'] != 0]
        
        # Store the processed dataframe
        enzyme_info_dfs[flux_col] = condition_df
        
        logger.debug(f"  {flux_col}: {len(condition_df)} rows after filtering")
    
    logger.info(f"Created enzyme info dataframes for {len(enzyme_info_dfs)} conditions")
    return enzyme_info_dfs


# Constants
_SECONDS_PER_HOUR = 3600.0
_FLUX_H_TO_UMOL_MIN = 1000.0 / 60.0  # SA unit conversion factor: mmol/gDCW/h -> umol/gDCW/min

# Reaction classes handled by the grouped (mass-fraction) calculator
GROUPED_CLASSES = ('isoenzyme', 'complex', 'mixed')


def calculate_kapp_homomeric(enzyme_protein_info_dfs: dict):
    """
    Calculate kapp for homomeric enzymes for each condition.

    For a homomeric enzyme the apparent turnover number is the classic ratio of
    flux to molar enzyme abundance:

        kcat_app [s^-1] = (v / 3600) / E

    where ``v`` is flux [mmol/gDCW/h] and ``E`` is ``protein_mmol_gdcw``
    [mmol/gDCW]. The mass fraction (``MW * E``) and specific activity are also
    computed so that the homomeric output shares one schema with the grouped
    (isoenzyme/complex) output.

    Parameters:
        enzyme_protein_info_dfs: dict
            Dictionary with structure {condition: dataframe}
    
    Returns:
        dict: Same structure as input but with added kapp columns in each dataframe
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
        df_copy = df_copy.drop_duplicates(subset=["gene", "SMILES", "rxn"])

        # Convert negative fluxes to positive (take absolute value)
        df_copy['flux_value'] = df_copy['flux_value'].abs()
        
        # COBRA fluxes are in mmol/gDW*h, convert to mmol/gDW*s
        df_copy['flux_value_per_sec'] = df_copy['flux_value'] / _SECONDS_PER_HOUR  # mmol/gDW*s
        
        # Calculate kcat_app: flux (mmol/gDW*s) / enzyme concentration (mmol/gDCW) = kcat (1/s)
        # Handle division by zero by replacing with NaN
        df_copy['kcat_app'] = df_copy['flux_value_per_sec'] / df_copy['protein_mmol_gdcw']
        
        # Replace infinite values with NaN
        df_copy['kcat_app'] = df_copy['kcat_app'].replace([float('inf'), float('-inf')], float('nan'))

        # Mass fraction [mg/gDCW] = MW [g/mol = mg/mmol] * abundance [mmol/gDCW]
        df_copy['mass_fraction'] = df_copy['molecular_weight'] * df_copy['protein_mmol_gdcw']
        df_copy['MW_total'] = df_copy['molecular_weight']
        # Specific activity [umol/mg/min] = flux [umol/gDCW/min] / mass_fraction [mg/gDCW]
        df_copy['SA_app'] = (df_copy['flux_value'] * _FLUX_H_TO_UMOL_MIN) / df_copy['mass_fraction']
        df_copy['SA_app'] = df_copy['SA_app'].replace([float('inf'), float('-inf')], float('nan'))

        # Unified identity / metadata columns (single-gene enzyme)
        df_copy['enzyme_key'] = df_copy['sequence']
        df_copy['genes'] = df_copy['gene']
        df_copy['sequences'] = df_copy['sequence']
        df_copy['n_contributing'] = 1
        if 'enzyme_class' not in df_copy.columns:
            df_copy['enzyme_class'] = 'homomeric'
        
        # Count valid kcat_app values
        valid_kcat = df_copy['kcat_app'].notna().sum()
        logger.debug(f"  {condition_name}: {initial_rows} → {len(df_copy)} rows (homomeric+dedup), {valid_kcat} valid kcat_app")
        
        # Store the processed dataframe
        kapp_results[condition_name] = df_copy
    
    logger.info("Completed kcat_app calculation for all conditions")
    return kapp_results


def _select_contributing_genes(enzyme_class, dnf, genes_with_data):
    """
    Decide which genes contribute to the integrated-enzyme mass fraction,
    applying the missing-data policy.

    Parameters
    ----------
    enzyme_class : str
        'isoenzyme', 'complex', or 'mixed'.
    dnf : list[frozenset[str]]
        DNF of the reaction GPR (OR of AND-groups).
    genes_with_data : set[str]
        Genes that have both a molecular weight and a measured abundance.

    Returns
    -------
    set[str]
        Genes to include in the sums. Empty -> the enzyme should be dropped.

    Policy
    ------
    - complex (single AND-group): require every subunit to have data; otherwise
      drop (a complex needs all subunits).
    - isoenzyme (OR of single genes): include every alternative that has data
      (an unmeasured/unexpressed isozyme contributes nothing).
    - mixed (OR of AND-groups): keep only OR-branches whose subunits all have
      data, and include the union of their genes.
    """
    if not dnf:
        return set()

    if enzyme_class == 'isoenzyme':
        return {g for branch in dnf for g in branch if g in genes_with_data}

    if enzyme_class == 'complex':
        branch = set().union(*dnf)
        return set(branch) if branch.issubset(genes_with_data) else set()

    if enzyme_class == 'mixed':
        kept = set()
        for branch in dnf:
            if branch.issubset(genes_with_data):
                kept |= set(branch)
        return kept

    # Fallback: anything with data
    return set(genes_with_data)


def calculate_kapp_grouped(enzyme_protein_info_dfs: dict, enzyme_classes=GROUPED_CLASSES):
    """
    Calculate kapp for isoenzymes and complexes using the Davidi et al. (2016)
    specific-activity method.

    For each reaction the full GPR is treated as one integrated enzyme. The
    enzyme mass fraction sums the contribution of every polypeptide chain:

        mass_fraction [mg/gDCW] = sum_g ( MW_g [mg/mmol] * abundance_g [mmol/gDCW] )

    The apparent specific activity and the (active-sites = 1, stoichiometry = 1)
    turnover number are then:

        SA_app   [umol/mg/min] = (v * 1000 / 60) / mass_fraction
        kcat_app [s^-1]        = (v / 3600) * MW_total / mass_fraction
                                = SA_app * MW_total / 60

    where ``v`` is flux [mmol/gDCW/h] and ``MW_total = sum_g MW_g`` is the total
    holoenzyme weight. This reduces exactly to the homomeric formula when a
    single gene is involved.

    Parameters
    ----------
    enzyme_protein_info_dfs : dict
        Dictionary {condition: dataframe}; each frame must carry
        'enzyme_id', 'enzyme_class', 'GPR', 'gene', 'molecular_weight',
        'protein_mmol_gdcw', 'flux_value', 'SMILES', 'Direction'.
    enzyme_classes : iterable[str]
        Subset of {'isoenzyme', 'complex', 'mixed'} to compute.

    Returns
    -------
    dict
        {condition: dataframe} with one row per integrated enzyme-substrate
        entry and the columns added by this function.
    """
    from src.enzyme_classifier import gpr_to_dnf

    target_classes = [c for c in enzyme_classes if c in GROUPED_CLASSES]
    kapp_results = {}

    if not target_classes:
        logger.info("No grouped enzyme classes requested; skipping grouped kapp.")
        return {name: None for name in enzyme_protein_info_dfs}

    group_keys = ['enzyme_id', 'SMILES', 'Direction']

    for condition_name, enzyme_df in enzyme_protein_info_dfs.items():
        if enzyme_df is None:
            kapp_results[condition_name] = None
            continue

        df = enzyme_df.copy()
        df = df[df['enzyme_class'].isin(target_classes)]
        if len(df) == 0:
            kapp_results[condition_name] = df
            logger.debug(f"  {condition_name}: no grouped enzymes")
            continue

        df['flux_value'] = df['flux_value'].abs()

        # Cache the DNF for each distinct GPR rule encountered in this condition
        dnf_cache = {}

        records = []
        for _, group in df.groupby(group_keys, dropna=False):
            first = group.iloc[0]
            gpr_rule = first['GPR']
            enzyme_class = first['enzyme_class']

            if gpr_rule not in dnf_cache:
                dnf_cache[gpr_rule] = gpr_to_dnf(gpr_rule)
            dnf = dnf_cache[gpr_rule]

            # One (MW, abundance) per distinct gene in this group
            gene_rows = group.drop_duplicates(subset=['gene'])
            gene_data = {}
            for _, row in gene_rows.iterrows():
                mw = row['molecular_weight']
                ab = row['protein_mmol_gdcw']
                if pd.notna(mw) and pd.notna(ab):
                    gene_data[row['gene']] = (float(mw), float(ab))

            contributing = _select_contributing_genes(enzyme_class, dnf, set(gene_data))
            if not contributing:
                continue

            mass_fraction = sum(gene_data[g][0] * gene_data[g][1] for g in contributing)
            mw_total = sum(gene_data[g][0] for g in contributing)
            if not (mass_fraction > 0):
                continue

            flux_value = float(first['flux_value'])
            flux_per_sec = flux_value / _SECONDS_PER_HOUR
            sa_app = (flux_value * _FLUX_H_TO_UMOL_MIN) / mass_fraction
            kcat_app = flux_per_sec * mw_total / mass_fraction

            sorted_genes = sorted(contributing)
            seq_map = (
                gene_rows.set_index('gene')['sequence'].to_dict()
                if 'sequence' in gene_rows.columns else {}
            )
            sequences = ';'.join(str(seq_map.get(g, '')) for g in sorted_genes)

            records.append({
                'enzyme_key': first['enzyme_id'],
                'enzyme_id': first['enzyme_id'],
                'enzyme_class': enzyme_class,
                'gpr_class': first.get('gpr_class'),
                'genes': ','.join(sorted_genes),
                'sequences': sequences,
                'sequence': sequences,
                'gene': ','.join(sorted_genes),
                'rxn': first['rxn'],
                'subsystem': first.get('subsystem'),
                'SMILES': first['SMILES'],
                'Direction': first['Direction'],
                'flux_value': flux_value,
                'flux_value_per_sec': flux_per_sec,
                'FVA_lower': first.get('FVA_lower'),
                'FVA_upper': first.get('FVA_upper'),
                'mass_fraction': mass_fraction,
                'MW_total': mw_total,
                'protein_mmol_gdcw': mass_fraction / mw_total if mw_total else float('nan'),
                'n_components': int(first.get('n_components', len(gene_data))),
                'n_contributing': len(contributing),
                'SA_app': sa_app,
                'kcat_app': kcat_app,
            })

        grouped_df = pd.DataFrame(records)
        kapp_results[condition_name] = grouped_df
        logger.debug(
            f"  {condition_name}: {len(grouped_df)} integrated enzymes "
            f"(classes: {target_classes})"
        )

    logger.info("Completed grouped (isoenzyme/complex) kcat_app calculation for all conditions")
    return kapp_results


def calculate_kapp(enzyme_protein_info_dfs: dict, enzyme_classes=('homomeric',)):
    """
    Calculate kapp for the requested enzyme classes and combine
    the per-condition results into a single dict of dataframes.

    Parameters
    ----------
    enzyme_protein_info_dfs : dict
        {condition: dataframe} after proteomics mapping.
    enzyme_classes : iterable[str]
        Any of {'homomeric', 'isoenzyme', 'complex', 'mixed'}.

    Returns
    -------
    dict
        {condition: dataframe} with kapp columns, concatenating the homomeric
        and grouped results that were requested.
    """
    classes = list(enzyme_classes)
    logger.info(f"Calculating kapp for enzyme classes: {classes}")

    homomeric_results = (
        calculate_kapp_homomeric(enzyme_protein_info_dfs)
        if 'homomeric' in classes else {}
    )

    grouped_classes = [c for c in classes if c in GROUPED_CLASSES]
    grouped_results = (
        calculate_kapp_grouped(enzyme_protein_info_dfs, enzyme_classes=grouped_classes)
        if grouped_classes else {}
    )

    combined = {}
    for condition_name in enzyme_protein_info_dfs:
        frames = []
        homo = homomeric_results.get(condition_name)
        grp = grouped_results.get(condition_name)
        if homo is not None and len(homo) > 0:
            frames.append(homo)
        if grp is not None and len(grp) > 0:
            frames.append(grp)
        if frames:
            combined[condition_name] = pd.concat(frames, ignore_index=True)
        else:
            combined[condition_name] = None

    return combined


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


# The kcat_app threshold filter is class-agnostic.
evaluate_kapp = evaluate_kapp_homomeric


def get_kapp_dataframe(kapp_filtered_results: dict) -> pd.DataFrame:
    """
    Collate kapp results from all conditions into a single wide-format DataFrame.

    Each unique (sequence, SMILES, gene, rxn, subsystem) combination becomes one row.
    Columns are the kcat_app value for each condition

    Parameters
    ----------
    kapp_filtered_results : dict
        Dictionary with structure {condition_name: pd.DataFrame} as returned by
        ``evaluate_kapp_homomeric``.

    Returns
    -------
    pd.DataFrame
        Wide-format DataFrame with one row per enzyme-substrate-reaction entry and
        one column per growth condition.
    """
    index_cols = ['enzyme_key', 'enzyme_class', 'sequence', 'SMILES', 'gene', 'rxn', 'subsystem']
    long_frames = []

    for condition_name, df in kapp_filtered_results.items():
        if df is None or len(df) == 0:
            logger.debug(f"  Skipping empty/None dataframe for condition: {condition_name}")
            continue

        available_index = [c for c in index_cols if c in df.columns]
        if 'kcat_app' not in df.columns:
            logger.warning(f"  No 'kcat_app' for condition '{condition_name}', skipping")
            continue

        subset = df[available_index + ['kcat_app']].copy()
        subset['_condition'] = condition_name
        long_frames.append(subset)

    if not long_frames:
        logger.warning("No valid data found across all conditions; returning empty DataFrame")
        return pd.DataFrame(columns=index_cols)

    long_df = pd.concat(long_frames, ignore_index=True)
    logger.debug(f"Combined long DataFrame: {len(long_df)} rows across {len(long_frames)} conditions")

    available_index = [c for c in index_cols if c in long_df.columns]

    # Pivot: one column per condition, values are kcat_app
    kapp_wide = long_df.pivot_table(
        index=available_index,
        columns='_condition',
        values='kcat_app',
        aggfunc='first',
    ).reset_index()

    # Remove the column-level name left by pivot_table
    kapp_wide.columns.name = None

    n_pairs = len(kapp_wide)
    n_conditions = len(kapp_filtered_results)
    logger.info(
        f"Collation summary: {n_pairs} unique enzyme-substrate entries "
        f"across {n_conditions} conditions"
    ) 

    return kapp_wide


def get_kmax_homomeric(kapp_results: dict):
    """
    Get the maximum kapp of each enzyme-substrate pair across all conditions.

    The enzyme is identified by ``enzyme_key`` (the protein sequence for
    homomeric enzymes, the integrated-enzyme id for isoenzymes/complexes), so
    this handles all enzyme classes uniformly. ``kcat_app_max`` is the maximum
    turnover (s^-1) and ``SA_app_max`` is the maximum specific activity
    (umol/mg/min) across conditions.

    Parameters:
        kapp_results: dict
            Dictionary with structure {condition: dataframe}
    Returns:
        kmax_results: pd.DataFrame
            DataFrame with the maximum kcat_app value for each enzyme-substrate
            pair plus metadata and SA_app_max.
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
        return pd.DataFrame(columns=['enzyme_key', 'sequence', 'SMILES', 'kcat_app_max', 'condition_max'])
    
    # Concatenate all dataframes
    combined_df = pd.concat(all_dataframes, ignore_index=True)
    logger.debug(f"Combined dataframe has {len(combined_df)} total entries")

    group_cols = ['enzyme_key', 'SMILES']

    # Representative row = the condition with the maximum turnover number
    kmax_results = (
        combined_df.loc[combined_df.groupby(group_cols)['kcat_app'].idxmax()]
        .reset_index(drop=True)
    )

    # Maximum specific activity across conditions (may differ from the kcat-max
    # condition when the contributing subunit set varies between conditions)
    sa_max = (
        combined_df.groupby(group_cols)['SA_app'].max().reset_index()
        .rename(columns={'SA_app': 'SA_app_max'})
        if 'SA_app' in combined_df.columns else None
    )

    # Select and rename relevant columns for output
    output_columns = [
        'enzyme_key', 'enzyme_class', 'sequence', 'sequences', 'genes',
        'SMILES', 'Direction', 'kcat_app', 'SA_app', 'source_condition',
        'gene', 'rxn', 'flux_value', 'FVA_upper', 'FVA_lower',
        'protein_mmol_gdcw', 'mass_fraction', 'MW_total',
        'n_components', 'n_contributing', 'subsystem',
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

    # Attach the across-condition SA maximum
    if sa_max is not None:
        kmax_results = kmax_results.merge(sa_max, on=group_cols, how='left')
    
    # Sort by kcat_app_max in descending order
    kmax_results = kmax_results.sort_values('kcat_app_max', ascending=False).reset_index(drop=True)
    
    logger.info(
        f"Found kmax for {len(kmax_results)} enzyme-substrate pairs "
        f"(range: {kmax_results['kcat_app_max'].min():.2e} to {kmax_results['kcat_app_max'].max():.2e} s⁻¹)"
    )
    
    return kmax_results


# Class-agnostic alias
get_kmax = get_kmax_homomeric


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
            kmax_results[['enzyme_key', 'SMILES', 'kcat_app_max']],
            on=['enzyme_key', 'SMILES'],
            how='left'
        )
        
        # Calculate eta = kcat_app / kcat_app_max
        df_with_eta['eta'] = df_with_eta['kcat_app'] / df_with_eta['kcat_app_max']
        
        # Replace infinite values with NaN
        df_with_eta['eta'] = df_with_eta['eta'].replace([float('inf'), float('-inf')], float('nan'))
        
        # Store the dataframe with eta
        kapp_results_with_eta[condition_name] = df_with_eta
        
        # Collect eta values for variance calculation
        eta_data = df_with_eta[['enzyme_key', 'SMILES', 'eta']].copy()
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
    variance_metrics = all_eta_df.groupby(['enzyme_key', 'SMILES'])['eta'].agg([
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
        variance_metrics[['enzyme_key', 'SMILES', 'eta_mean', 'eta_stdev', 'eta_min', 'eta_max', 'eta_cv']],
        on=['enzyme_key', 'SMILES'],
        how='left'
    )
    
    # Log summary statistics
    logger.info(
        f"Eta calculation complete: {len(kmax_with_variance)} pairs, "
        f"mean η={variance_metrics['eta_mean'].mean():.3f}, CV={variance_metrics['eta_cv'].mean():.3f}"
    )
    
    return kapp_results_with_eta, kmax_with_variance
