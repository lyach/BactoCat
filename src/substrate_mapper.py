'''
TO DO:
1. Make a function that will map the enzyme substrates in a GEM. Use the dataset_kinGEMs.py script as a reference.

- INPUT: GEM model
- OUTPUT: DataFrame with the substrates for each reaction in the GEM, if its reversible or not, and the SMILES string of the substrate.
'''

from misc.dataset_kinGEMs import map_metabolites
from misc.dataset_kinGEMs import get_substrate_metabolites

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
