import argparse
import csv
from pathlib import Path
import cobra 
import pandas as pd
import numpy as np
from loguru import logger

# ---------------------------------------------------------------------
# 1. SUPPLEMENTARY MAPPING CSV (The "Database")
# ---------------------------------------------------------------------

def create_supplementary_csv(model: cobra.Model, output_dir: Path) -> Path:
    """
    Create AIDA compound → GEM exchange reaction mapping CSV.
    This structure supports 1-to-many mappings (salts) and stoichiometry.
    """
    # Mapping: { AIDA_Name: { GEM_Rxn: Multiplier } }
    mapping_rules = {
        # --- Single Matches ---

        # --- Carbon Sources ---
        "Glucose (mM)": {"EX_glc__D_e": 1.0},

        # --- Amino Acids ---
        "Alanine (mM)": {"EX_ala__L_e": 1.0},
        "Arginine/HCl (mM)": {"EX_arg__L_e": 1.0},
        "Asparagine/H2O (mM)": {"EX_asn__L_e": 1.0},
        "AsparticAcid (mM)": {"EX_asp__L_e": 1.0},
        "GlutamicAcid/HCl (mM)": {"EX_glu__L_e": 1.0},
        "Glutamine (mM)": {"EX_gln__L_e": 1.0},
        "Glycine (mM)": {"EX_gly_e": 1.0},
        "Histidine/HCl/H2O (mM)": {"EX_his__L_e": 1.0},
        "Isoleucine (mM)": {"EX_ile__L_e": 1.0},
        "Leucine (mM)": {"EX_leu__L_e": 1.0},
        "Lysine/HCl (mM)": {"EX_lys__L_e": 1.0},
        "Methionine (mM)": {"EX_met__L_e": 1.0},
        "Phenylalanine (mM)": {"EX_phe__L_e": 1.0},
        "Proline (mM)": {"EX_pro__L_e": 1.0},
        "Serine (mM)": {"EX_ser__L_e": 1.0},
        "Threonine (mM)": {"EX_thr__L_e": 1.0},
        "Tryptophan (mM)": {"EX_trp__L_e": 1.0},
        "Tyrosine (mM)": {"EX_tyr__L_e": 1.0},
        "Valine (mM)": {"EX_val__L_e": 1.0},
        "Thiamine/HCl (mM)": {"EX_thm_e": 1.0},
        "Pyridoxine (mM)": {"EX_pydxn_e": 1.0},
        # or "Cystine/HCl/H2O (mM)": {"EX_cys__L_e": 2.0},  

        # --- Ambiguous / Salt Matches ---

        "Cystine/HCl/H2O (mM)": {"EX_cys__L_e": 2.0, "EX_cl_e": 2.0},  # Cystine is 2x Cysteine + 2x Cl

        # --- Nitrogen, Phosphorus, Sulfur ---
        "(NH4)2SO4 (mM)": {"EX_nh4_e": 2.0, "EX_so4_e": 1.0},
        "NH4Cl (mM)": {"EX_nh4_e": 1.0, "EX_cl_e": 1.0},
        
        "K2HPO4 (mM)": {"EX_pi_e": 1.0, "EX_k_e": 2.0},
        "KH2PO4 (mM)": {"EX_pi_e": 1.0, "EX_k_e": 1.0},
        "Na2HPO4 (mM)": {"EX_pi_e": 1.0, "EX_na1_e": 2.0},

        # --- Magnesium and Calcium ---
        "MgCl2/6H2O (mM)": {"EX_mg2_e": 1.0, "EX_cl_e": 2.0},
        "MgSO4/7H2O (mM)": {"EX_mg2_e": 1.0, "EX_so4_e": 1.0},
        "CaCl2/2H2O (mM)": {"EX_ca2_e": 1.0, "EX_cl_e": 2.0},
        "CaSO4/2H2O (mM)": {"EX_ca2_e": 1.0, "EX_so4_e": 1.0},

        # --- Iron and Trace Metals ---
        "FeSO4/7H2O (mM)": {"EX_fe2_e": 1.0, "EX_so4_e": 1.0},
        "ZuSO4/7H2O (mM)": {"EX_zn2_e": 1.0, "EX_so4_e": 1.0}, # Zinc (typo in AIDA)
        "CuSO4/5H2O (mM)": {"EX_cu2_e": 1.0, "EX_so4_e": 1.0},
        "Na2MoO4/2H2O (mM)": {"EX_mobd_e": 1.0, "EX_na1_e": 2.0},
        
        # --- Others ---
        "KCl (mM)": {"EX_k_e": 1.0, "EX_cl_e": 1.0},
        "NaCl (mM)": {"EX_na1_e": 1.0, "EX_cl_e": 1.0},
        "Na3C6H5O7/2H2O (mM)": {"EX_cit_e": 1.0, "EX_na1_e": 3.0}, # Citrate
        "Na2S2O3/5H2O (mM)": {"EX_tsul_e": 1.0, "EX_na1_e": 2.0},  # Thiosulfate
    }

    logger.info(f"Total AIDA compounds in mapping: {len(mapping_rules)}")
    valid_rxns = {r.id for r in model.reactions}
    out_csv = output_dir / "ecoli_aida_direct_mappings.csv"

    with out_csv.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["AIDA_compound", "gem_rxn", "stoichiometry", "match_type"])

        for comp, rxns in mapping_rules.items():
            # Classify based on number of mapped reactions
            match_type = "single_match" if len(rxns) == 1 else "ambiguous"
            
            for rxn, coeff in rxns.items():
                if rxn in valid_rxns:
                    writer.writerow([comp, rxn, coeff, match_type])
                else:
                    writer.writerow([comp, rxn, coeff, "missing_match"])
                    logger.warning(f"Reaction {rxn} not found in GEM for {comp}")

    return out_csv

# ---------------------------------------------------------------------
# 2. MATRIX BUILDER (New Helper)
# ---------------------------------------------------------------------

def build_mapping_matrix(mapping_csv: Path, media_columns: list) -> pd.DataFrame:
    """Constructs M matrix: Rows=AIDA, Cols=GEM reactions."""
    map_df = pd.read_csv(mapping_csv)
    
    # Identify which AIDA compounds from our mapping actually exist in the CSV
    present_aida_cols = [c for c in map_df["AIDA_compound"].unique() if c in media_columns]
    unique_gem_rxns = map_df["gem_rxn"].unique().tolist()

    # Initialize zero matrix
    M = pd.DataFrame(0.0, index=present_aida_cols, columns=unique_gem_rxns)

    for _, row in map_df.iterrows():
        if row["AIDA_compound"] in M.index:
            M.loc[row["AIDA_compound"], row["gem_rxn"]] = row["stoichiometry"]
    
    return M

# ---------------------------------------------------------------------
# 3. MEDIA + GROWTH CSV (Transformation Engine)
# ---------------------------------------------------------------------

def create_media_growth_csv(aida_dir: Path, mapping_csv: Path, output_dir: Path) -> Path:
    media_csv = aida_dir / "medium_composition.csv"
    growth_csv = aida_dir / "growth_data.csv"

    media_df = pd.read_csv(media_csv)
    growth_df = pd.read_csv(growth_csv)

    # Preserve original
    orig_path = aida_dir / "media_composition_original.csv"
    media_df.to_csv(orig_path, index=False)
    logger.info(f"Preserved original media CSV → {orig_path}")

    # Normalize media column names
    media_df.columns = (
        media_df.columns
        .str.replace("\n", " ", regex=False)
        .str.replace(r"\s+", " ", regex=True)
        .str.strip()
    )

    # ------------------------------------------------------------
    # Build and Apply Mapping Matrix
    # ------------------------------------------------------------
    M = build_mapping_matrix(mapping_csv, media_df.columns.tolist())
    
    # Identify what was NOT renamed/mapped
    not_renamed_cols = [
        col for col in media_df.columns 
        if col not in M.index and col != "Condition ID"
    ]
    
    logger.info(f"Mapped {len(M.index)} AIDA columns into {len(M.columns)} GEM reactions.")
    logger.info(f"{len(not_renamed_cols)} media columns were not mapped.")
    logger.info(f"Columns not renamed: {not_renamed_cols}")

    # Use Matrix Multiplication to transform the data (Handles Summing & Stoichiometry)
    A = media_df[M.index].fillna(0.0)
    gem_values_df = A.dot(M)
    
    # Combine back with Condition ID
    media_gem_df = pd.concat([media_df[["Condition ID"]], gem_values_df], axis=1)
    gem_rxns = M.columns.tolist()

    # ------------------------------------------------------------
    # Deduplicate by unique media composition
    # ------------------------------------------------------------
    logger.info(f"Media compositions before dropping duplicates → {media_gem_df.shape[0]}")
    media_gem_df = media_gem_df.drop_duplicates(subset=gem_rxns).reset_index(drop=True)
    logger.info(f"Unique media compositions after dropping duplicates → {media_gem_df.shape[0]}")

    # ------------------------------------------------------------
    # Get average growth data (r_info filtering)
    # ------------------------------------------------------------
    logger.info(f"Total growth entries in CSV: {len(growth_df)}")
    logger.info(f"Entries with r_info = 1: {len(growth_df[growth_df['r_info'] == 1])}")
    logger.info(f"Entries with r_info = 2: {len(growth_df[growth_df['r_info'] == 2])}")

    # Only keep r_info = 0
    growth_df = growth_df[growth_df["r_info"] == 0]
    logger.info(f"Growth entries remaining (r_info = 0): {len(growth_df)}")

    avg_growth = (
        growth_df
        .groupby("Condition ID", as_index=False)
        .agg(avg_growth=("r", "mean"))
    )

    # ------------------------------------------------------------
    # Merge media + growth
    # ------------------------------------------------------------
    final_df = pd.merge(media_gem_df, avg_growth, on="Condition ID", how="inner")

    # Reorder columns
    final_df = final_df[["Condition ID", "avg_growth"] + gem_rxns]

    out_csv = output_dir / "ecoli_aida_media_growth.csv"
    final_df.to_csv(out_csv, index=False)

    logger.info(f"Wrote final media + growth CSV → {out_csv}")
    logger.info(f"Final unique conditions → {final_df.shape[0]}")

    return out_csv

# ---------------------------------------------------------------------
# 4. MAIN (Includes License-Fix)
# ---------------------------------------------------------------------

def main():
    # Force GLPK to bypass expired Gurobi license before loading model
    cobra.Configuration().solver = "glpk"

    parser = argparse.ArgumentParser()
    parser.add_argument("--model_path", required=True, type=Path)
    parser.add_argument("--aida_dir", required=True, type=Path)
    parser.add_argument("--output_dir", required=True, type=Path)
    args = parser.parse_args()

    logger.info(f"Loading GEM → {args.model_path}")
    model = cobra.io.read_sbml_model(str(args.model_path))

    mapping_csv = create_supplementary_csv(model, args.output_dir)
    create_media_growth_csv(args.aida_dir, mapping_csv, args.output_dir)

if __name__ == "__main__":
    main()