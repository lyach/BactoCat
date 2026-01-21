"""
AIDA → E. coli GEM media & growth CSV builder

Outputs:
1) ecoli_aida_direct_mappings.csv
2) media_composition_original.csv
3) ecoli_aida_media_growth.csv
"""

import argparse
import csv
from pathlib import Path

import cobra 
import pandas as pd
from loguru import logger


# ---------------------------------------------------------------------
# SUPPLEMENTARY MAPPING CSV
# ---------------------------------------------------------------------

def create_supplementary_csv(model: cobra.Model, output_dir: Path) -> Path:
    """
    Create AIDA compound → GEM exchange reaction mapping CSV.
    """

    single_matches = {
        "Glucose (mM)": "EX_glc__D_e",

        "Alanine (mM)": "EX_ala__L_e",
        "Arginine/HCl (mM)": "EX_arg__L_e",
        "Asparagine/H2O (mM)": "EX_asn__L_e",
        "AsparticAcid (mM)": "EX_asp__L_e",
        "GlutamicAcid/HCl (mM)": "EX_glu__L_e",
        "Glutamine (mM)": "EX_gln__L_e",
        "Glycine (mM)": "EX_gly_e",
        "Histidine/HCl/H2O (mM)": "EX_his__L_e",
        "Isoleucine (mM)": "EX_ile__L_e",
        "Leucine (mM)": "EX_leu__L_e",
        "Lysine/HCl (mM)": "EX_lys__L_e",
        "Methionine (mM)": "EX_met__L_e",
        "Phenylalanine (mM)": "EX_phe__L_e",
        "Proline (mM)": "EX_pro__L_e",
        "Serine (mM)": "EX_ser__L_e",
        "Threonine (mM)": "EX_thr__L_e",
        "Tryptophan (mM)": "EX_trp__L_e",
        "Tyrosine (mM)": "EX_tyr__L_e",
        "Valine (mM)": "EX_val__L_e",

        "Thiamine/HCl (mM)": "EX_thm_e",
        "Pyridoxine (mM)": "EX_pydxn_e",

        "(NH4)2SO4 (mM)": "EX_nh4_e",  # Maps concentration directly to Nitrogen
        "NH4Cl (mM)": "EX_nh4_e",       # Maps concentration directly to Nitrogen
        "K2HPO4 (mM)": "EX_pi_e",       # Maps to Phosphorus
        "MgSO4/7H2O (mM)": "EX_mg2_e",  # Maps to Magnesium

    }

    # log the number of single matches
    logger.info(f"Length of single matches: {len(single_matches)}")

    ambiguous_matches = {
        #"K2HPO4 (mM)": "EX_pi_e",
        "KH2PO4 (mM)": "EX_k_e",
        "Na2HPO4 (mM)": "EX_na1_e",

        #"(NH4)2SO4 (mM)": ["EX_nh4_e", "EX_so4_e"], 
        #"NH4Cl (mM)": ["EX_nh4_e", "EX_cl_e"],

        "Cystine/HCl/H2O (mM)": ["EX_cys__L_e"],

        "KCl (mM)": ["EX_k_e"],
        "NaCl (mM)": ["EX_na1_e"],

        "MgCl2/6H2O (mM)": ["EX_mg2_e"],
        #"MgSO4/7H2O (mM)": ["EX_mg2_e", "EX_so4_e"],

        "CaSO4/2H2O (mM)": ["EX_ca2_e", "EX_so4_e"],
        "CaCl2/2H2O (mM)": "EX_ca2_e",

        "FeSO4/7H2O (mM)": ["FE2abcpp", "EX_fe2_e"],

        "ZnSO4/7H2O (mM)": "EX_zn2_e",
        
        "Na3C6H5O7/2H2O (mM)": "EX_cit_e",
        "Na2S2O3/5H2O (mM)": "EX_tsul_e",
        "CuSO4/5H2O (mM)": "EX_cu2_e",
        "Na2MoO4/2H2O (mM)": "EX_mobd_e",

    }

    # log the number of ambiguous matches
    logger.info(f"Length of ambiguous matches: {len(ambiguous_matches)}")

    # log the number of single and ambiguous matches combined
    total_matches = len(single_matches) + len(ambiguous_matches)
    logger.info(f"Total number of AIDA single + ambiguous matches: {total_matches}")

    valid_rxns = {r.id for r in model.reactions}
    output_dir.mkdir(parents=True, exist_ok=True)
    out_csv = output_dir / "ecoli_aida_direct_mappings.csv"

    with out_csv.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["AIDA_compound", "gem_rxn", "match_type"])

        for comp, rxn in single_matches.items():
            if rxn in valid_rxns:
                writer.writerow([comp, rxn, "single_match"])
            else:
                writer.writerow([comp, "", "missing_match"])

        for comp, rxns in ambiguous_matches.items():
            for rxn in rxns:
                if rxn in valid_rxns:
                    writer.writerow([comp, rxn, "ambiguous"])

    logger.info(f"Wrote supplementary mapping CSV → {out_csv}")
    return out_csv

# ---------------------------------------------------------------------
# Build Mapping Matrix for salts
# ---------------------------------------------------------------------


def build_mapping_matrix(mapping_csv: Path, media_columns: list) -> pd.DataFrame:
    """
    Converts the mapping CSV into a Stoichiometric Matrix (M).
    Rows: AIDA Compounds
    Cols: GEM Reactions
    """
    map_df = pd.read_csv(mapping_csv)
    
    # Only include AIDA compounds that actually appear in our media CSV
    present_aida_cols = [c for c in map_df["AIDA_compound"].unique() if c in media_columns]
    unique_gem_rxns = map_df["gem_rxn"].unique()

    # Initialize the Matrix M with zeros
    M = pd.DataFrame(0.0, index=present_aida_cols, columns=unique_gem_rxns)

    # Populate M with stoichiometry
    for _, row in map_df.iterrows():
        comp = row["AIDA_compound"]
        rxn = row["gem_rxn"]
        val = row["stoichiometry"]
        if comp in M.index:
            M.loc[comp, rxn] = val
            
    logger.info(f"Built mapping matrix: {M.shape[0]} AIDA compounds -> {M.shape[1]} GEM reactions")
    return M

# ---------------------------------------------------------------------
# MEDIA + GROWTH CSV
# ---------------------------------------------------------------------

def create_media_growth_csv(
    aida_dir: Path,
    mapping_csv: Path,
    output_dir: Path,
) -> Path:
    """
    Final CSV:
    [Condition ID, avg_growth, EX_glc__D_e, EX_ala__L_e, ...]
    One row per unique media composition.
    """

    media_csv = aida_dir / "medium_composition.csv"
    growth_csv = aida_dir / "growth_data.csv"

    media_df = pd.read_csv(media_csv)
    growth_df = pd.read_csv(growth_csv)

    # Preserve original
    orig_path = aida_dir / "media_composition_original.csv"
    media_df.to_csv(orig_path, index=False)
    logger.info(f"Preserved original media CSV → {orig_path}")

    # ------------------------------------------------------------
    # Normalize media column names 
    # ------------------------------------------------------------
    media_df.columns = (
        media_df.columns
        .str.replace("\n", " ", regex=False)
        .str.replace(r"\s+", " ", regex=True)
        .str.strip()
    )

    # ------------------------------------------------------------
    # Load mapping
    # ------------------------------------------------------------
    mapping_df = pd.read_csv(mapping_csv)
    mapping_df = mapping_df[mapping_df["match_type"] == "single_match"]

    comp_to_rxn = dict(
        zip(mapping_df["AIDA_compound"], mapping_df["gem_rxn"])
    )

    # ------------------------------------------------------------
    # Rename media columns → GEM rxns
    # ------------------------------------------------------------
    rename_cols = {
        col: comp_to_rxn[col]
        for col in media_df.columns
        if col in comp_to_rxn
    }

    # CALCULATE NOT RENAMED BEFORE FILTERING
    # We find columns that are NOT in our rename dict and NOT the ID column
    not_renamed_cols = [
        col for col in media_df.columns 
        if col not in rename_cols and col != "Condition ID"
    ]

    # Log the counts and the specific names
    num_renamed = len(rename_cols)
    num_not_renamed = len(not_renamed_cols)
    
    logger.info(f"Renamed {num_renamed} media columns to GEM exchange reactions.")
    logger.info(f"{num_not_renamed} media columns were not renamed.")
    logger.info(f"Columns not renamed: {not_renamed_cols}")

    # NOW rename and filter the dataframe
    media_df = media_df.rename(columns=rename_cols)
    gem_rxns = list(rename_cols.values())
    media_df = media_df[["Condition ID"] + gem_rxns]

    # ------------------------------------------------------------
    # Deduplicate by unique media composition
    # ------------------------------------------------------------
    
    logger.info(f"Number of media compositions before dropping duplicates→ {media_df.shape[0]}")

    media_df = (
        media_df
        .drop_duplicates(subset=gem_rxns) # remove duplicate media compositions
        .reset_index(drop=True)
    )
    logger.info(f"Unique media compositions after dropping duplicates → {media_df.shape[0]}")

    # ------------------------------------------------------------
    # Get average growth data
    # ------------------------------------------------------------

    # 1. Log the absolute total first
    logger.info(f"Total growth entries in CSV: {len(growth_df)}")

    # 2. Log the breakdown (for your information)
    logger.info(f"Entries with r_info = 1: {len(growth_df[growth_df['r_info'] == 1])}")
    logger.info(f"Entries with r_info = 2: {len(growth_df[growth_df['r_info'] == 2])}")

    # 3. Perform the filter
    growth_df = growth_df[growth_df["r_info"].isin([0])] # only keep r_info = 0

    avg_growth = (
        growth_df
        .groupby("Condition ID", as_index=False)
        .agg(avg_growth=("r", "mean")) # average growth rate (r)
    )

    # 4. Log the result after filtering
    logger.info(f"Growth entries remaining (r_info = 0): {len(growth_df)}")
    # ------------------------------------------------------------
    # Merge media + growth
    # ------------------------------------------------------------
    final_df = pd.merge(
        media_df,
        avg_growth,
        on="Condition ID",
        how="inner" # only keep conditions with growth data (r_info = 0)
    )

    # Reorder columns
    final_df = final_df[
        ["Condition ID", "avg_growth"] + gem_rxns
    ]

    out_csv = output_dir / "ecoli_aida_media_growth.csv"
    final_df.to_csv(out_csv, index=False)

    logger.info(f"Wrote final media + growth CSV → {out_csv}")
    logger.info(f"Final unique conditions → {final_df.shape[0]}")

    return out_csv


# ---------------------------------------------------------------------
# MAIN
# ---------------------------------------------------------------------

def main():
    cobra.Configuration().solver = "glpk"

    parser = argparse.ArgumentParser()
    parser.add_argument("--model_path", required=True, type=Path)
    parser.add_argument("--aida_dir", required=True, type=Path)
    parser.add_argument("--output_dir", required=True, type=Path)
    args = parser.parse_args()

    logger.info(f"Loading GEM → {args.model_path}")
    model = cobra.io.read_sbml_model(args.model_path)

    mapping_csv = create_supplementary_csv(model, args.output_dir)
    create_media_growth_csv(args.aida_dir, mapping_csv, args.output_dir)


if __name__ == "__main__":
    main()
