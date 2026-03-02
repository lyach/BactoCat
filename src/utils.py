"""
Module Description: utils.py

Purpose: 
Utility functions for the BactoCat pipeline.
"""

import pandas as pd
import warnings
from rdkit import Chem
warnings.filterwarnings('ignore', category=RuntimeWarning)


# =============================================================================
# Modeling functions
# =============================================================================

def

# =============================================================================
# Functions for Aida dataset preprocessing
# =============================================================================

def prepare_aida_dataset(df_pred_path: pd.DataFrame, 
                         df_cond_path: pd.DataFrame) -> pd.DataFrame:
    """
    Prepare the Aida dataset to be used in the BactoCat pipeline.
    """
    # Load
    df_pred = pd.read_csv(df_pred_path)
    print(f"Loaded Aida predictions df with shape {df_pred.shape}")
    
    df_cond = pd.read_csv(df_cond_path)
    print(f"Loaded Aida conditions df with shape {df_cond.shape}")
    
    # TO DO: Check that the number of conditions in df_cond 
    # is the same as the number of conditions in df_pred
    # and populate the condition_id column in df_pred with the Condition ID column from df_cond
    # which has the og condition number from Aida & Ying
    
    # TEMPORARY: Add placeholder id
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
# Functions for specific kcat datasets
# =============================================================================

def load_kcat_dataset_ecoli(CPIPred_dir, CatPred_dir, EnzyExtract_dir) -> pd.DataFrame:
    """
    Load in vitro kcat datasets for E. coli.
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
# Others
# =============================================================================

def canonicalize(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        return Chem.MolToSmiles(mol, canonical=True)
    except:
        return None