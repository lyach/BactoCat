"""
paxdb_mapper.py

Purpose: 
Map PaxDB protein abundances to an enzyme dataframe by gene identifier and calculate protein concentrations.

Overview: 
This module provides functions to:
1. Extract gene IDs from PaxDB `string_external_id` (e.g., "511145.b0440" -> "b0440").
2. Aggregate abundances per gene (handles duplicates).
3. Map abundances to the enzyme dataframe as `protein_ppm`.
4. Calculate molecular weights from protein sequences using BioPython.
5. Convert protein abundances to molar concentrations (mmol/gDCW).
"""

import pandas as pd
from Bio.SeqUtils.ProtParam import ProteinAnalysis


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
    
    # Calculate protein_mol_gdcw: protein_ppm * p_total / molecular_weight
    # where p_total is the total protein content in g/gDCW
    # molecular_weight is in g/mol
    enz_mapped["protein_mol_gdcw"] = (
        enz_mapped["protein_fraction"] * p_total / enz_mapped["molecular_weight"]
    )
    
    # mol/gDCW to mmol/gDCW
    enz_mapped["protein_mmol_gdcw"] = enz_mapped["protein_mol_gdcw"] * 1000

    return enz_mapped
