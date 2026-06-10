"""
utils.py

Utility functions for the BactoCat pipeline.
"""

import re
import numpy as np
import pandas as pd
import warnings
import cobra
import os

from scipy import stats

from loguru import logger
from cobra import flux_analysis
from rdkit import Chem

warnings.filterwarnings('ignore', category=RuntimeWarning)


# =============================================================================
# Constants
# =============================================================================


M9_MEDIA = {'EX_glc__D_e', 'EX_so4_e', 'EX_o2_e', 'EX_nh4_e', 'EX_pi_e', 'EX_h2o_e',
    'EX_h_e', 'EX_k_e', 'EX_mg2_e', 'EX_na1_e', 'EX_ca2_e', 'EX_cl_e', 'EX_cu2_e', 'EX_fe2_e', 
    'EX_fe3_e', 'EX_mn2_e', 'EX_mobd_e', 'EX_ni2_e', 'EX_zn2_e', 'EX_cobalt2_e'}

M9_MEDIA_BOUND = 1000.0

CARBON_SOURCE_MAP = {
    'glucose':     'EX_glc__D_e',
    'galactose':   'EX_gal_e',
    'acetate':     'EX_ac_e',
    'pyruvate':    'EX_pyr_e',
    'fumarate':    'EX_fum_e',
    'succinate':   'EX_succ_e',
    'glucosamine': 'EX_gam_e',
    'glycerol':    'EX_glyc_e',
    'mannose':     'EX_man__D_e',
    'xylose':      'EX_xyl__D_e',
    'fructose':    'EX_fru_e',
}


# =============================================================================
# Modeling functions
# =============================================================================


def get_constrained_growth(model: cobra.Model, 
                           medium_df: pd.DataFrame,
                           biomass_rxn: str = "BIOMASS_Ec_iML1515_core_75p37M",
                           method: str = "pFBA") -> pd.DataFrame:
    """
    Get the constrained growth rate of a model given a medium dataframe.
    
    Parameters
    ----------
    model : cobra.Model
        The COBRA model to optimize
    medium_df : pd.DataFrame
        The medium dataframe to use for the optimization
    biomass_rxn : str, optional
        The biomass reaction to use for the optimization
    method : str, optional
        The method to use for the optimization (FBA or pFBA)

    Returns
    -------
    pd.DataFrame
        A dataframe with columns 'condition_id', 'growth_exp', and 'growth_sim'
    """
    non_medium_cols = {'condition_id', 'avg_growth'}
    results = []

    for _, row in medium_df.iterrows():
        condition_id = row['condition_id']
        avg_growth = row['avg_growth']

        medium_dict = {
            col: row[col]
            for col in medium_df.columns
            if col not in non_medium_cols
        }

        with model:
            from src.kapp_builder import modify_reaction_bounds

            modify_reaction_bounds(model, medium_dict, medium_upper_bound=False, verbose=False)

            if method == "FBA":
                solution = model.optimize()
            elif method == "pFBA":
                solution = flux_analysis.pfba(model)
            else:
                raise ValueError(f"Invalid method '{method}'. Must be 'FBA' or 'pFBA'.")

            growth_sim = solution.fluxes[biomass_rxn]

        results.append({
            'condition_id': condition_id,
            'growth_exp': avg_growth,
            'growth_sim': growth_sim
        })

    return pd.DataFrame(results)


# =============================================================================
# Functions for AMN dataset preprocessing
# =============================================================================


def prepare_amn_dataset(df_pred_path: pd.DataFrame) -> pd.DataFrame:
    """
    Prepare the AMN dataset to be used in the BactoCat pipeline.
    """
    # Load
    df_pred = pd.read_csv(df_pred_path)
    print(f"Loaded AMN predictions df with shape {df_pred.shape}")
    
    # Add condition id
    df_pred['condition_id'] = [f'cond{i+1}' for i in range(len(df_pred))]
    
    # Rename GR_AVG column to avg_growth
    df_pred.rename(columns={'GR_AVG': 'avg_growth'}, inplace=True)
    
    # Rename all columns (except for condition_id and avg_growth) to remove the suffix '_i'
    protected_cols = ['condition_id', 'avg_growth']
    
    new_columns = []
    for col in df_pred.columns:
        if col not in protected_cols and col.endswith('_i'):
            new_columns.append(col[:-2]) # Remove last 2 chars
        else:
            new_columns.append(col)
    
    df_pred.columns = new_columns
    
    # Move condition_id and avg_growth columns to the front
    cols = list(df_pred.columns)
    for col in ['avg_growth', 'condition_id']:
        if col in cols:
            cols.insert(0, cols.pop(cols.index(col)))
    
    # Reorder columns
    df_pred = df_pred[cols]
    
    return df_pred


# =============================================================================
# Functions for Davidi dataset preprocessing
# =============================================================================


def davidi_condition_to_id(condition_name: str) -> str:
    """
    Normalise a Davidi-style growth condition IDs.
    """
    s = condition_name.lower()
    s = re.sub(r'[+@=]', '_', s)
    s = re.sub(r'(?<=\d)\.(?=\d)', '', s)
    s = re.sub(r'[^a-z0-9]+', '_', s)
    s = s.strip('_')
    return s

def prepare_davidi_dataset(growth_conditions_path: str) -> pd.DataFrame:
    """
    Convert Davidi growth conditions CSV into the BactoCat conditions format.

    All M9_MEDIA exchange reactions are set to M9_MEDIA_BOUND (1000.0).
    For non-glucose carbon sources, EX_glc__D_e is set to 0 and the
    appropriate exchange reaction is added at M9_MEDIA_BOUND.
    """
    df = pd.read_csv(growth_conditions_path)
    logger.info(f"Loaded Davidi growth conditions with shape {df.shape}")

    rows = []
    for _, row in df.iterrows():
        growth_cond = str(row['growth condition'])
        growth_rate = float(row['growth rate [h-1]'])
        carbon_source = str(row['media (M9 plus)']).strip().lower()

        condition_id = davidi_condition_to_id(growth_cond)

        new_row = {'condition_id': condition_id, 'avg_growth': growth_rate}

        for ex in M9_MEDIA:
            new_row[ex] = M9_MEDIA_BOUND

        carbon_source_ex = CARBON_SOURCE_MAP.get(carbon_source)
        if carbon_source_ex is None:
            logger.warning(f"Unknown carbon source '{carbon_source}' for condition '{growth_cond}'; EX_glc__D_e kept open.")
        elif carbon_source_ex != 'EX_glc__D_e':
            new_row['EX_glc__D_e'] = 0.0
            new_row[carbon_source_ex] = M9_MEDIA_BOUND

        rows.append(new_row)

    result_df = pd.DataFrame(rows)

    fixed_cols = ['condition_id', 'avg_growth']
    media_cols = sorted(c for c in result_df.columns if c not in fixed_cols)
    result_df = result_df[fixed_cols + media_cols]

    logger.info(f"Prepared Davidi conditions df with shape {result_df.shape}")
    
    return result_df


# =============================================================================
# Functions for specific kcat datasets
# =============================================================================


def load_kcat_dataset_ecoli(CPIPred_dir: str, 
                            CatPred_dir: str, 
                            EnzyExtract_dir: str) -> pd.DataFrame:
    """
    Load in vitro kcat datasets for E. coli.
    
    Parameters:
        CPIPred_dir: Path to the CPIPred CSV file
        CatPred_dir: Path to the CatPred CSV file
        EnzyExtract_dir: Path to the EnzyExtract parquet file
    
    Returns:
        pd.DataFrame: DataFrame containing the loaded kcat datasets
    """
    # Load
    CPIPred_df = pd.read_csv(CPIPred_dir)
    CatPred_df = pd.read_csv(CatPred_dir)
    EnzyExtract_df = pd.read_parquet(EnzyExtract_dir)
    
    # Keep only E coli data
    CPIPred_df = CPIPred_df[CPIPred_df['organism'].str.contains("Escherichia coli", case=False, na=False)]
    CatPred_df = CatPred_df[(CatPred_df['taxonomy_id'] == 562) | (CatPred_df['taxonomy_id'] == 83333)]
    EnzyExtract_df = EnzyExtract_df[EnzyExtract_df['organism'].str.contains("Escherichia coli", case=False, na=False)]
    
    # Keep and rename useful columns
    CPIPred_df = CPIPred_df[["SEQ", "CMPD_SMILES", "kcat"]]
    CPIPred_df = CPIPred_df[CPIPred_df['kcat'].notna()]
    CPIPred_df.rename(columns={"SEQ": "sequence", "CMPD_SMILES": "SMILES", "kcat": "kcat_CPIPred"}, inplace=True)

    CatPred_df = CatPred_df[["sequence", "reactant_smiles", "value"]]
    CatPred_df = CatPred_df[CatPred_df['value'].notna()]
    CatPred_df.rename(columns={'reactant_smiles': 'SMILES', "value": "kcat_CatPred"}, inplace=True)
    
    EnzyExtract_df = EnzyExtract_df[["sequence", "smiles", "kcat_value"]]
    EnzyExtract_df = EnzyExtract_df.dropna(subset=['kcat_value', 'sequence', 'smiles'])
    EnzyExtract_df.rename(columns={"kcat_value": "kcat_EnzyExtract", "smiles": "SMILES"}, inplace=True)


    #df_kcat = pd.concat([CPIPred_df, CatPred_df])
        
    return CPIPred_df, CatPred_df, EnzyExtract_df


def process_catpred_smiles(df: pd.DataFrame, smiles_col: str = 'reactant_smiles') -> pd.DataFrame:
    """
    Process CatPred SMILES entries to extract main substrates, filtering out cofactors.
    
    Parameters
    ----------
    df : pd.DataFrame
        Input dataframe containing SMILES entries
    smiles_col : str
        Column name containing the SMILES entries to process
        
    Returns
    -------
    pd.DataFrame
        New dataframe with additional 'SMILES' column containing individual 
        substrate SMILES, with one row per substrate. Original rows with multiple 
        substrates are expanded into multiple rows.
    """
    
    # Define common cofactors to filter out
    cofactors = {
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
    }
    
    # Create a copy of the dataframe
    result_df = df.copy()
    
    # List to store processed rows
    processed_rows = []
    
    for idx, row in result_df.iterrows():
        smiles_entry = row[smiles_col]
        
        # Split SMILES by '.' to get individual molecules
        molecules = [mol.strip() for mol in str(smiles_entry).split('.')]
        
        # Filter out cofactors and keep only substrates
        substrates = []
        for mol in molecules:
            if mol and mol not in cofactors:
                # Additional filtering for very simple molecules that are likely cofactors
                # Skip molecules that are too simple (less than 3 non-hydrogen atoms)
                if _is_likely_substrate(mol):
                    substrates.append(mol)
        
        # If no substrates found, keep the original entry
        if not substrates:
            row_copy = row.copy()
            row_copy['SMILES'] = smiles_entry
            processed_rows.append(row_copy)
        else:
            # Create a row for each substrate
            for substrate in substrates:
                row_copy = row.copy()
                row_copy['SMILES'] = substrate
                processed_rows.append(row_copy)
    
    # Convert back to DataFrame
    result_df = pd.DataFrame(processed_rows)
    
    # Drop duplicates
    result_df = result_df.drop_duplicates(subset=["sequence", "SMILES"])

    result_df.reset_index(drop=True, inplace=True)
    
    return result_df


def _is_likely_substrate(smiles: str) -> bool:
    """
    Determine if a SMILES string represents a likely substrate (not a cofactor).
    
    This function uses heuristics to distinguish between substrates and cofactors
    based on molecular complexity and common cofactor patterns.
    
    Parameters
    ----------
    smiles : str
        SMILES string representing a molecule
        
    Returns
    -------
    bool
        True if the molecule is likely a substrate, False if likely a cofactor
    """
    
    # Additional known cofactors not caught by simple string matching
    simple_cofactors = {
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
    
    if smiles in simple_cofactors:
        return False
    
    # Count non-hydrogen atoms (rough complexity measure)
    # This is a simplified approach - in reality you'd use rdkit for proper atom counting
    non_h_chars = sum(1 for c in smiles if c.isupper() and c not in ['H'])
    
    # If very few heavy atoms, likely a cofactor
    if non_h_chars < 3:
        return False
    
    # Additional heuristics can be added here
    # For now, if it passes the above filters, consider it a substrate
    return True


# =============================================================================
# Metrics
# =============================================================================

def _log_transform(x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    mask = (x > 0) & (y > 0)
    return np.log10(x[mask]), np.log10(y[mask])


def _r2_score(x: np.ndarray, y: np.ndarray, log: bool = False) -> tuple[float, float, int]:
    if log:
        x, y = _log_transform(x, y)
    # R² from linear regression
    slope, intercept, r_value, _, _ = stats.linregress(x, y)
    r2 = float(r_value ** 2)
    return r2, float('nan'), len(x)


def _pearson_r(x: np.ndarray, y: np.ndarray, log: bool = False) -> tuple[float, float, int]:
    if log:
        x, y = _log_transform(x, y)
    r, pval = stats.pearsonr(x, y)
    return float(r), float(pval), len(x)


def _spearmans_rho(x: np.ndarray, y: np.ndarray, log: bool = False) -> tuple[float, float, int]:
    rho, pval = stats.spearmanr(x, y)
    return float(rho), float(pval), len(x)


def _kendall_tau(x: np.ndarray, y: np.ndarray, log: bool = False) -> tuple[float, float, int]:
    tau, pval = stats.kendalltau(x, y)
    return float(tau), float(pval), len(x)


def _rmsd(x: np.ndarray, y: np.ndarray, log: bool = False) -> tuple[float, float, int]:
    if log:
        x, y = _log_transform(x, y)
    rmsd = float(np.sqrt(np.mean((x - y) ** 2)))
    return rmsd, float('nan'), len(x)


def get_metrics(df: pd.DataFrame, value1: str, value2: str, log: bool = True) -> pd.DataFrame:
    """
    Compute comparison metrics between two columns.

    Parameters
    ----------
    df : pd.DataFrame
    value1 : str
        Column name for the first set of values.
    value2 : str
        Column name for the second set of values.
    log : bool
        If True, passes log=True to each metric that supports it (R², Pearson r,
        RMSD). Spearman and Kendall are unaffected by log scaling.

    Returns
    -------
    pd.DataFrame with columns: metric, value, p_value, samples
    """
    sub = df[[value1, value2]].dropna()
    x = sub[value1].to_numpy(dtype=float)
    y = sub[value2].to_numpy(dtype=float)

    _metrics = [
        ('R²',            _r2_score),
        ('Pearson r',     _pearson_r),
        ("Spearman's ρ",  _spearmans_rho),
        ("Kendall's τ",   _kendall_tau),
        ('RMSD',          _rmsd),
    ]

    rows = []
    for name, func in _metrics:
        stat, pval, n = func(x, y, log=log)
        rows.append({'metric': name, 'value': stat, 'p_value': pval, 'samples': n})

    return pd.DataFrame(rows)


# =============================================================================
# ETA analysis
# =============================================================================


def preprocess_in_vitro_dataset(path: str, dataset_name: str) -> pd.DataFrame:
    """Standardizes columns, types, and string formatting."""
    df = pd.read_parquet(path)

    df.columns = (
        df.columns
        .str.strip()
        .str.lower()
    )
      
    if dataset_name.upper() == "ENZYEXTRACT":
        # Adjust these keys based on your specific parquet's headers
        rename_map = {
            'sequence': 'sequence', 
            'smiles': 'SMILES',
            'kcat_value': 'kcat_in_vitro'
        }
        df = df.rename(columns=rename_map)

        # Drop rows where essential merge keys are actually missing (NaN)
        df = df.dropna(subset=['sequence', 'SMILES', 'kcat_in_vitro'])
        
        return df[['sequence', 'SMILES', 'kcat_in_vitro']]
    
    elif dataset_name.upper() == "DAVIDI":
        # TO DO
        pass
    
    return df

def get_eta_in_vitro(in_vitro_kcat_path: str, kmax_results: pd.DataFrame, kmax_path: str, dataset_name: str):
    """
    Calculate eta (k_in_vitro / kmax) for a specific dataset.
    Based on the original get_eta but modified for single-source in vitro comparison.

    Paremeters:
        in_vitro_kcat_path: Path to the in vitro kcat dataset (parquet)
        kmax_results: DataFrame containing kmax results with 'sequence', 'SMILES',
                        'kcat_app_max', and 'subsystem' columns 
        kmax_path: Path to the original kmax results (used for saving output in same directory)
        dataset_name: Name of the in vitro dataset (used for logging and output file naming)
   
     Returns:
        pd.DataFrame: DataFrame containing the calculated eta values and related information.
    """
    logger.info(f"Calculating eta using {dataset_name} data...")

    # 1. Preprocess the external dataset
    in_vitro_df = preprocess_in_vitro_dataset(in_vitro_kcat_path, dataset_name)

    # Re-drop any that failed canonicalization
    kmax_results = kmax_results.dropna(subset=['SMILES'])
    in_vitro_df = in_vitro_df.dropna(subset=['SMILES'])
    kmax_results['SMILES'] = kmax_results['SMILES'].apply(canonicalize)
    in_vitro_df['SMILES'] = in_vitro_df['SMILES'].apply(canonicalize)

    # 2. Merge with kmax_results
    merged_df = pd.merge(
        kmax_results[['sequence', 'SMILES', 'kcat_app_max', 'subsystem']],
        in_vitro_df,
        on=['sequence', 'SMILES'],
        how='inner' 
    )

    logger.info(f"Number of matches between kmax and {dataset_name}: {len(merged_df)}")

    # 3. Calculate eta = kcat_app_max / kcat_in_vitro
    merged_df['eta'] = merged_df['kcat_app_max'] / merged_df['kcat_in_vitro']
    merged_df['eta'] = merged_df['eta'].replace([float('inf'), float('-inf')], float('nan'))

    # 4. Save to CSV in the same folder as kmax_path
    output_dir = os.path.dirname(kmax_path)
    output_file = os.path.join(output_dir, f'eta_{dataset_name}.csv')

    # Final column ordering
    final_cols = ['sequence', 'SMILES', 'kcat_app_max', 'kcat_in_vitro', 'subsystem', 'eta']
    output_df = merged_df[final_cols]
        
    output_df.to_csv(output_file, index=False)
    logger.info(f"Analysis saved to: {output_file}")

    return output_df


# =============================================================================
# Others
# =============================================================================


def canonicalize(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        return Chem.MolToSmiles(mol, canonical=True)
    except:
        return None
    
    
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