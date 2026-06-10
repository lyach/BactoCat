"""
proteomics_mapper.py

Module to map protein abundances to an enzyme dataframe by gene 
identifierand calculate protein concentrations in mmol/gDCW.
"""

import pandas as pd
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from loguru import logger

from src.utils import davidi_condition_to_id


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
    ----------
    paxdb_df : pd.DataFrame
        PaxDB dataframe with 'string_external_id' and 'abundance' columns
    df_enzymes : pd.DataFrame
        Enzyme dataframe with 'gene' and 'sequence' columns
    p_total : float
        Total protein content in g/gDCW
    
    Returns
    -------
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

    # If multiple entries per gene, take the mean
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
    
    # Calculate protein_mol_gdcw: protein_ppm * p_total / molecular_weight
    # where p_total is the total protein content in g/gDCW
    # molecular_weight is in g/mol
    enz_mapped["protein_mol_gdcw"] = (
        enz_mapped["protein_fraction"] * p_total / enz_mapped["molecular_weight"]
    )
    
    # mol/gDCW to mmol/gDCW
    enz_mapped["protein_mmol_gdcw"] = enz_mapped["protein_mol_gdcw"] * 1000

    return enz_mapped


def paxdb_protein_mapping(enzyme_info_dfs: dict, paxdb_path: str, p_total: float):
    """
    Apply PaxDB protein mapping across all enzyme info dataframes.

    Parameters
    ----------
    enzyme_info_dfs : dict
        Dictionary with condition names as keys and enzyme dataframes as values
    paxdb_path : str
        Path to PaxDB TSV file
    p_total : float
        Total protein content value (g protein / g DCW)

    Returns
    -------
    dict
        Dictionary structure: {condition_name: mapped_dataframe}
    """
    logger.info(f"Loading PaxDB data from: {paxdb_path}")
    paxdb_df = pd.read_csv(
        paxdb_path,
        sep="\t",
        comment="#",
        header=None,
        names=["gene_name", "string_external_id", "abundance"],
    )
    logger.debug(f"PaxDB data loaded: {len(paxdb_df)} rows")

    enzyme_protein_info_dfs = {}

    logger.info(f"Processing {len(enzyme_info_dfs)} conditions with p_total={p_total}")

    for condition_name, enzyme_df in enzyme_info_dfs.items():
        try:
            mapped_df = map_paxdb_to_gene(
                paxdb_df=paxdb_df,
                df_enzymes=enzyme_df,
                p_total=p_total,
            )

            enzyme_protein_info_dfs[condition_name] = mapped_df
            logger.debug(
                f"  {condition_name}: {len(mapped_df)} rows, "
                f"{mapped_df['protein_ppm'].notna().sum()} with protein data"
            )

        except Exception as e:
            logger.error(f"Error processing {condition_name}: {str(e)}")
            enzyme_protein_info_dfs[condition_name] = None

    return enzyme_protein_info_dfs


def map_davidi_to_gene(proteomics_df: pd.DataFrame, df_enzymes: pd.DataFrame, condition_col: str) -> pd.DataFrame:
    """
    Map Davidi proteomics abundances to enzymes by gene ID for a specific condition.

    Parameters
    ----------
    proteomics_df : pd.DataFrame
        Proteomics dataframe with 'gene' column and condition columns
    df_enzymes : pd.DataFrame
        Enzyme dataframe with a 'gene' column
    condition_col : str
        Column name in proteomics_df corresponding to the target condition

    Returns
    -------
    pd.DataFrame
        Copy of df_enzymes with a new 'protein_mmol_gdcw' column (NaN if no match)
    """
    if "gene" not in df_enzymes.columns:
        raise KeyError("df_enzymes missing required column: 'gene'")
    if condition_col not in proteomics_df.columns:
        raise KeyError(f"Condition column '{condition_col}' not found in proteomics data")

    prot = proteomics_df[["gene", condition_col]].copy()
    prot["protein_mmol_gdcw"] = pd.to_numeric(prot[condition_col], errors="coerce")
    prot = prot[["gene", "protein_mmol_gdcw"]]

    return df_enzymes.copy().merge(prot, on="gene", how="left")


def specific_proteome_mapping(enzyme_info_dfs: dict, proteomics_path: str) -> dict:
    """
    Apply condition-specific proteome mapping across all enzyme info dataframes.

    Parameters
    ----------
    enzyme_info_dfs : dict
        Dictionary with condition names as keys and enzyme dataframes as values.
        Keys must already be normalised via `davidi_condition_to_id`.
    proteomics_path : str
        Path to proteomics CSV with columns:
        'bnumber', <gene_name>, <condition_1>, <condition_2>, ...
        where condition columns contain protein abundance in mmol/gDCW.

    Returns
    -------
    dict
        Dictionary structure: {condition_name: mapped_dataframe}
    """
    logger.info(f"Loading proteomics data from: {proteomics_path}")
    prot_df = pd.read_csv(proteomics_path)

    # Standardise gene identifier column
    prot_df = prot_df.rename(columns={prot_df.columns[0]: "gene"})
    # Drop gene name column (second column)
    prot_df = prot_df.drop(columns=[prot_df.columns[1]])
    # Normalise condition column names
    condition_cols = list(prot_df.columns[1:])
    prot_df = prot_df.rename(columns={col: davidi_condition_to_id(col) for col in condition_cols})

    logger.debug(f"Proteomics data loaded: {len(prot_df)} genes, {len(prot_df.columns) - 1} conditions")

    enzyme_protein_info_dfs = {}

    logger.info(f"Processing {len(enzyme_info_dfs)} conditions")

    for condition_name, enzyme_df in enzyme_info_dfs.items():
        try:
            mapped_df = map_davidi_to_gene(
                proteomics_df=prot_df,
                df_enzymes=enzyme_df,
                condition_col=condition_name,
            )
            enzyme_protein_info_dfs[condition_name] = mapped_df
            logger.debug(
                f"  {condition_name}: {len(mapped_df)} rows, "
                f"{mapped_df['protein_mmol_gdcw'].notna().sum()} with protein data"
            )
        except Exception as e:
            logger.error(f"Error processing {condition_name}: {str(e)}")
            enzyme_protein_info_dfs[condition_name] = None

    return enzyme_protein_info_dfs