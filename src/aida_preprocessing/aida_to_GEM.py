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
        "Aspartic Acid (mM)": "EX_asp__L_e",
        "Glutamic Acid/HCl (mM)": "EX_glu__L_e",
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


    }

    print(f"Length of single matches: {len(single_matches)}")

    ambiguous_matches = {
        "K2HPO4 (mM)": "EX_pi_e",
        "KH2PO4 (mM)": "EX_k_e",
        "Na2HPO4 (mM)": "EX_na1_e",

        "(NH4)2SO4 (mM)": ["EX_nh4_e", "EX_so4_e"], # check 
        "NH4Cl (mM)": ["EX_nh4_e", "EX_cl_e"],

        "Cystine/HCl/H2O (mM)": ["CYSTL", "cysi__L","cys__L"],

        "KCl (mM)": ["Kabcpp", "Kt2pp", "EX_k_e"],
        "NaCl (mM)": ["NAt3_1p5pp", "NAt3_2pp", "EX_na1_e"],

        "MgCl2/6H2O (mM)": ["MG2uabcpp", "EX_mg2_e"],
        "MgSO4/7H2O (mM)": ["EX_mg2_e", "EX_so4_e"],

        "CaSO4/2H2O (mM)": ["EX_ca2_e", "EX_so4_e"],
        "CaCl2/2H2O (mM)": "EX_ca2_e",

        "FeSO4/7H2O (mM)": ["FE2abcpp", "EX_fe2_e"],

        "Folic Acid (mM)": ["FTHFD", "DHFR"],
        "Aminobenzoic Acid (mM)": ["ADCL", "DHPS2"],

        "H3BO3 (mM)": ["hbo3_e"],

        "Riboflavin (mM)": ["ribflv_c"],

        "ZnSO4/7H2O (mM)": "EX_zn2_e",
        
        "Na3C6H5O7/2H2O (mM)": "EX_cit_e",
        "Na2S2O3/5H2O (mM)": "EX_tsul_e",
        "CuSO4/5H2O (mM)": "EX_cu2_e",
        "Na2MoO4/2H2O (mM)": "EX_mobd_e",

    }

    print(f"Length of ambiguous matches: {len(ambiguous_matches)}")

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

    media_df = media_df.rename(columns=rename_cols)

    gem_rxns = list(rename_cols.values())

    media_df = media_df[["Condition ID"] + gem_rxns]

    # ------------------------------------------------------------
    # Deduplicate by unique media composition
    # ------------------------------------------------------------
    media_df = (
        media_df
        .drop_duplicates(subset=gem_rxns) # remove duplicate media compositions
        .reset_index(drop=True)
    )

    # ------------------------------------------------------------
    # Get average growth data
    # ------------------------------------------------------------
    growth_df = growth_df[growth_df["r_info"].isin([0])] # only keep r_info = 0

    avg_growth = (
        growth_df
        .groupby("Condition ID", as_index=False)
        .agg(avg_growth=("r", "mean")) # average growth rate (r)
    )

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
