"""
substrate_mapper.py

Purpose: 
Map substrates to SMILES structures using external databases.

Overview: 
The following functions in this file are attributed to the kinGEMs
project (https://github.com/ranaabarghout/kinGEMs/): 
- load_model  
- get_substrate_metabolites
- map_metabolites
- get_SMILES_with_retries
- get_SMILES_from_cactus
- get_PubChem_SMILES
- clean_metabolite_names
"""

import logging
import os
import random
import re
import time
import urllib
import urllib.error
from urllib.parse import quote
from urllib.request import urlopen
from tqdm import tqdm

from bioservices import UniProt
import cobra
from cobra.core import Reaction
from cobra.util.solver import set_objective
import numpy as np
import pandas as pd
import pubchempy as pcp

from .config import (
    CHEBI_COMPOUNDS,
    CHEBI_INCHI,
    METANETX_COMPOUNDS,
    METANETX_XREF,
    SEED_COMPOUNDS,
    TAXONOMY_IDS,
    BiGG_MAPPING,
    ensure_dir_exists,
)

def load_model(model_path):
    """
    Load a COBRA model from file.
    
    Parameters
    ----------
    model_path : str
        Path to the model file (SBML format)
        
    Returns
    -------
    cobra.Model
        The loaded model
    """
    try:
        model, errors = cobra.io.validate_sbml_model(model_path)
        if errors:
            print(f"Warning: Model has {len(errors)} validation errors")
        return cobra.io.read_sbml_model(model_path)
    except Exception as e:
        raise ValueError(f"Error loading model: {e}")

def get_substrate_metabolites(reaction):
    """
    Get the substrates (reactants) for a reaction.
    
    Parameters
    ----------
    reaction : cobra.Reaction
        The reaction to analyze
        
    Returns
    -------
    list
        List of substrate metabolite IDs
    """
    substrates = [met.id for met in reaction.reactants]
    return substrates

def map_metabolites(substrate_df, external_db_dir=None, max_retries=3, retry_delay=2):
    """
    Map metabolites to SMILES structures using external databases.
    
    Parameters
    ----------
    substrate_df : pandas.DataFrame
        DataFrame with substrate information
    external_db_dir : str, optional
        Directory containing external database files. If None, uses default.
    max_retries : int, optional
        Maximum number of retries for web service requests. Default is 3.
    retry_delay : int, optional
        Delay between retries in seconds. Default is 2.
        
    Returns
    -------
    pandas.DataFrame
        DataFrame with added SMILES information
    """
    import logging
    import os
    import re  # noqa: F811
    import time  # noqa: F401, F811
    from urllib.error import HTTPError  # noqa: F401

    import numpy as np
    import pandas as pd
    
    # Configure logging
    logging.basicConfig(level=logging.INFO, 
                        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    logger = logging.getLogger(__name__)
    
    # Load database files
    if external_db_dir is None:
        external_db_dir = os.path.dirname(BiGG_MAPPING)
    
    # Load database files
    BiGG_comps = pd.read_csv(BiGG_MAPPING)
    CHEBI_comps = pd.read_csv(CHEBI_COMPOUNDS, sep='\t')  # noqa: F841
    CHEBIInChI_comps = pd.read_csv(CHEBI_INCHI, sep='\t')
    MetaNetX_comps = pd.read_csv(METANETX_COMPOUNDS, sep='\t')
    MetaNetX_refcomps = pd.read_csv(METANETX_XREF, sep='\t')  # noqa: F841
    SEED_comps = pd.read_csv(SEED_COMPOUNDS, sep='\t')
    
    # Initialize new columns
    df = substrate_df.copy()
    df['SMILES'] = np.nan
    df['BiGG Name'] = np.nan
    df['DB Name'] = np.nan

    # Clean substrate names before processing
    df['Cleaned Substrate'] = df['Substrate partner'].apply(clean_metabolite_names)

    # Get unique substrates with their details
    unique_substrates_df = df[['Substrate partner', 'Cleaned Substrate']].drop_duplicates()
    logger.info(f"There are {len(unique_substrates_df)} substrates in the GEM.")

    # Prepare for tracking SMILES and names
    smiles_mapping = {}
    bigg_name_mapping = {}
    db_name_mapping = {}

    # Expand old_bigg_ids to create a list for comparison
    BiGG_comps['old_bigg_ids_list'] = BiGG_comps['old_bigg_ids'].str.split(';')
    BiGG_comps_expl = BiGG_comps.explode('old_bigg_ids_list')
    BiGG_comps_unique = BiGG_comps.reset_index(drop=True)
    BiGG_comps_expl_unique = BiGG_comps_expl.reset_index(drop=True)

    # Iterate through each unique substrate
    for _, row in tqdm(unique_substrates_df.iterrows(), total=len(unique_substrates_df), desc="Mapping substrates"):
        substrate = row['Substrate partner']
        cleaned_substrate = row['Cleaned Substrate']

        # Search in BiGG database
        bigg_hit = BiGG_comps_unique.loc[
            BiGG_comps_unique['bigg_id'].str.contains(fr'\b{re.escape(substrate)}\b', regex=True, case=False) |
            BiGG_comps_unique['universal_bigg_id'].str.contains(fr'\b{re.escape(substrate)}\b', regex=True, case=False) |
            BiGG_comps_expl_unique['old_bigg_ids_list'].str.contains(fr'\b{re.escape(substrate)}\b', regex=True, case=False)
        ].head(1)

        if not bigg_hit.empty:
            name = bigg_hit['name'].values[0]
            bigg_name_mapping[substrate] = name

            # Check MetaNetX for SMILES
            if not bigg_hit['MetaNetX'].isna().all():
                metanetx_hit = MetaNetX_comps[MetaNetX_comps['ID'] == bigg_hit['MetaNetX'].values[0]]
                if not metanetx_hit.empty:
                    smiles = metanetx_hit['SMILES'].values[0]
                    if pd.notna(smiles):
                        smiles_mapping[substrate] = smiles
                        db_name_mapping[substrate] = metanetx_hit['name'].values[0]
                        continue
            
            # Check SEED database
            if not bigg_hit['SEED'].isna().all():
                seed_hit = SEED_comps[SEED_comps['id'] == bigg_hit['SEED'].values[0]]
                if not seed_hit.empty:
                    smiles = seed_hit['smiles'].values[0] if 'smiles' in seed_hit else np.nan
                    name = seed_hit['name'].values[0]
                    if pd.notna(smiles):
                        smiles_mapping[substrate] = smiles
                        db_name_mapping[substrate] = name
                        continue
            
            # Check CHEBI database
            if not bigg_hit['CHEBI'].isna().all():
                try:
                    chebi_hit = CHEBIInChI_comps[CHEBIInChI_comps['CHEBI_ID'].str.contains(
                        fr'\b{re.escape(bigg_hit["CHEBI"].values[0])}\b', regex=True, case=False)]
                except AttributeError:
                    # If string methods fail, convert to string first
                    CHEBIInChI_comps['CHEBI_ID'] = CHEBIInChI_comps['CHEBI_ID'].astype(str)
                    chebi_hit = CHEBIInChI_comps[CHEBIInChI_comps['CHEBI_ID'].str.contains(
                        fr'\b{re.escape(bigg_hit["CHEBI"].values[0])}\b', regex=True, case=False)]
                
                if not chebi_hit.empty:
                    inchi = chebi_hit['InChI'].values[0]
                    inchi_hit = MetaNetX_comps[MetaNetX_comps['InChI'] == inchi]
                    if not inchi_hit.empty:
                        smiles = inchi_hit['SMILES'].values[0]
                        name = inchi_hit['name'].values[0]
                        if pd.notna(smiles):
                            smiles_mapping[substrate] = smiles
                            db_name_mapping[substrate] = name
                            continue
        
        # If no BiGG hit, check other databases directly
        else:            
            # Check MetaNetX directly
            metanetx_hit = MetaNetX_comps[MetaNetX_comps['ID'].str.contains(
                fr'\b{re.escape(substrate)}\b', regex=True, case=False)]
            if not metanetx_hit.empty:
                name = metanetx_hit['name'].values[0]
                smiles = metanetx_hit['SMILES'].values[0] if 'SMILES' in metanetx_hit.columns else np.nan
                if pd.notna(smiles):
                    smiles_mapping[substrate] = smiles
                    db_name_mapping[substrate] = name
                    continue
            
            # Check SEED
            seed_hit = SEED_comps[SEED_comps['id'].str.contains(
                fr'\b{re.escape(substrate)}\b', regex=True, case=False)]
            if not seed_hit.empty:
                name = seed_hit['name'].values[0]
                smiles = seed_hit['smiles'].values[0] if 'smiles' in seed_hit.columns else np.nan
                if pd.notna(smiles):
                    smiles_mapping[substrate] = smiles
                    db_name_mapping[substrate] = name
                    continue

    # If SMILES is still missing, try to get from web services
    missing_substrates = [sub for sub in df['Substrate partner'].unique() if sub not in smiles_mapping]
    

    for substrate in tqdm(missing_substrates, desc="Mapping missing substrates"):
        # First, check alternative databases directly
        metanetx_hit = MetaNetX_comps[MetaNetX_comps['ID'].str.contains(
            fr'\b{re.escape(substrate)}\b', regex=True, case=False)]
        if not metanetx_hit.empty:
            smiles = metanetx_hit['SMILES'].values[0] if 'SMILES' in metanetx_hit.columns else np.nan
            if pd.notna(smiles):
                smiles_mapping[substrate] = smiles
                db_name_mapping[substrate] = metanetx_hit['name'].values[0]
                continue

        seed_hit = SEED_comps[SEED_comps['id'].str.contains(
            fr'\b{re.escape(substrate)}\b', regex=True, case=False)]
        if not seed_hit.empty:
            smiles = seed_hit['smiles'].values[0] if 'smiles' in seed_hit.columns else np.nan
            if pd.notna(smiles):
                smiles_mapping[substrate] = smiles
                db_name_mapping[substrate] = seed_hit['name'].values[0]
                continue

        # Prepare for web service search
        cleaned_substrate = df.loc[df['Substrate partner'] == substrate, 'Cleaned Substrate'].iloc[0]
        
        # Try searches with both original and cleaned names
        names_to_try = [cleaned_substrate, substrate]
        
        for name in names_to_try:
            if pd.notna(name):
                # Try CIR with retries
                smiles = get_SMILES_with_retries(name, service='cactus', max_retries=max_retries, retry_delay=retry_delay)
                if smiles:
                    smiles_mapping[substrate] = smiles
                    break
                else:
                    # Try PubChem with retries
                    smiles_list = get_SMILES_with_retries(name, service='pubchem', max_retries=max_retries, retry_delay=retry_delay)
                    if smiles_list and len(smiles_list) > 0:
                        smiles_mapping[substrate] = smiles_list[0]
                        break

    # Apply mappings back to the dataframe
    for substrate, smiles in tqdm(smiles_mapping.items(), desc="Applying SMILES"):
        df.loc[df['Substrate partner'] == substrate, 'SMILES'] = smiles

    for substrate, bigg_name in tqdm(bigg_name_mapping.items(), desc="Applying BiGG names"):
        df.loc[df['Substrate partner'] == substrate, 'BiGG Name'] = bigg_name

    for substrate, db_name in tqdm(db_name_mapping.items(), desc="Applying DB names"):
        df.loc[df['Substrate partner'] == substrate, 'DB Name'] = db_name

    return df


def get_SMILES_with_retries(name, service='cactus', max_retries=3, retry_delay=2):
    """
    Get SMILES from various services with retry logic.
    
    Parameters
    ----------
    name : str
        Compound name to search
    service : str
        Service to use ('cactus' or 'pubchem')
    max_retries : int
        Maximum number of retries
    retry_delay : int
        Delay between retries in seconds
        
    Returns
    -------
    str or list
        SMILES string or list of SMILES strings
    """
    import logging
    import time  # noqa: F811
    
    logger = logging.getLogger(__name__)
    
    for attempt in range(1, max_retries + 1):
        try:
            if service == 'cactus':
                smiles = get_SMILES_from_cactus(name)
                return smiles
            elif service == 'pubchem':
                smiles_list = get_PubChem_SMILES(name)
                return smiles_list
        except Exception as e:
            if attempt < max_retries:
                time.sleep(retry_delay)
            else:
                return None
    
    return None


def get_SMILES_from_cactus(name):
    """
    Get SMILES from Chemical Identifier Resolver (CIR) with improved error handling.
    
    Parameters
    ----------
    name : str
        Compound name to search
        
    Returns
    -------
    str
        SMILES string
    """
    import logging
    import urllib.error
    import urllib.parse
    import urllib.request
    
    logger = logging.getLogger(__name__)
    
    # Properly encode the name for URL
    encoded_name = urllib.parse.quote(name)
    url = f"https://cactus.nci.nih.gov/chemical/structure/{encoded_name}/smiles"
    
    try:
        with urllib.request.urlopen(url, timeout=10) as response:
            smiles = response.read().decode('utf8')
            return smiles.strip()
    except urllib.error.HTTPError as e:
        return None
    except urllib.error.URLError as e:
        return None
    except Exception as e:
        return None


def get_PubChem_SMILES(name):
    """
    Get SMILES from PubChem with improved error handling.
    
    Parameters
    ----------
    name : str
        Compound name to search
        
    Returns
    -------
    list
        List of SMILES strings
    """
    import json
    import logging
    import time  # noqa: F811

    import requests
    
    logger = logging.getLogger(__name__)
    
    # First search for the compound to get CIDs
    search_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{urllib.parse.quote(name)}/cids/JSON"
    smiles_list = []
    
    try:
        search_response = requests.get(search_url, timeout=15)
        search_response.raise_for_status()
        search_data = search_response.json()
        
        if 'IdentifierList' in search_data and 'CID' in search_data['IdentifierList']:
            cids = search_data['IdentifierList']['CID']
            
            # Limit to first 3 results to avoid too many requests
            for cid in cids[:3]:
                # Small delay to avoid rate limiting
                time.sleep(0.5)
                
                # Get SMILES for each CID
                property_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/CanonicalSMILES/JSON"
                try:
                    prop_response = requests.get(property_url, timeout=15)
                    prop_response.raise_for_status()
                    prop_data = prop_response.json()
                    
                    if 'PropertyTable' in prop_data and 'Properties' in prop_data['PropertyTable']:
                        for compound in prop_data['PropertyTable']['Properties']:
                            if 'CanonicalSMILES' in compound:
                                smiles_list.append(compound['CanonicalSMILES'])
                except Exception as e:
                    continue
        
        return smiles_list
    except requests.exceptions.RequestException as e:
        return []
    except json.JSONDecodeError as e:
        return []
    except Exception as e:
        return []


def clean_metabolite_names(name):
    """
    Clean metabolite names for better matching in databases.
    
    Parameters
    ----------
    name : str
        Raw metabolite name
        
    Returns
    -------
    str
        Cleaned metabolite name
    """
    import re  # noqa: F811
    
    if pd.isna(name):
        return name
    
    # Convert to string if not already
    name = str(name)
    
    # Remove compartment suffix (e.g., _c, _e, _p)
    name = re.sub(r'_[a-z]$', '', name)
    
    # Remove common prefixes/suffixes
    name = re.sub(r'^(cpd|M_|m_)', '', name)
    
    # Replace underscores with spaces
    name = name.replace('_', ' ')
    
    # Remove concentration indicators like (e) or [e]
    name = re.sub(r'[\(\[].[^\)\]]*[\)\]]', '', name)
    
    # Remove charge indicators
    name = re.sub(r'[+-]\d*', '', name)
    
    # Strip whitespace
    name = name.strip()
    
    return name


def get_substrate_df(model, external_db_dir=None):
    """
    Generate a DataFrame mapping all substrates in a GEM.

    Parameters
    ----------
    model : cobra.Model
        The metabolic model (GEM)
    external_db_dir : str, optional
        Directory containing external database files for mapping metabolites

    Returns
    -------
    pandas.DataFrame
        DataFrame with columns:
        - 'Reaction': reaction ID
        - 'Substrate partner': substrate metabolite ID
        - 'Direction': 'forward' or 'reverse'
        - 'SMILES': SMILES string for the substrate
    """
    import pandas as pd

    rxn_data = []

    # Loop over all reactions
    for reaction in model.reactions:
        # Forward direction (reactants)
        substrates = get_substrate_metabolites(reaction)
        for substrate in substrates:
            rxn_data.append({
                "Reaction": reaction.id,
                "Substrate partner": substrate,
                "Direction": "forward"
            })

        # If reversible, add reverse direction (products as substrates)
        if getattr(reaction, "reversibility", False):
            products = [met.id for met in reaction.products]
            for product in products:
                rxn_data.append({
                    "Reaction": reaction.id,
                    "Substrate partner": product,
                    "Direction": "reverse"
                })

    # Create DataFrame 
    substrate_df = pd.DataFrame(rxn_data)

    # Map substrates to SMILES
    substrate_df_with_smiles = map_metabolites(substrate_df, external_db_dir)

    return substrate_df_with_smiles
