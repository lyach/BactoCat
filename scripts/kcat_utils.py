"""
Module Description: kcat_utils.py

Purpose: 
Utility functions for kcat visualization.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from typing import Tuple, Dict
import warnings
from rdkit import Chem
warnings.filterwarnings('ignore', category=RuntimeWarning)


# ============================================
# Plots for comparing kcat distributions
# ============================================

def compare_kcat_distribution(df1: pd.DataFrame, kcat_col1: str, 
                            df2: pd.DataFrame, kcat_col2: str,
                            label1: str = "Dataset 1", 
                            label2: str = "Dataset 2",
                            figsize: Tuple[int, int] = (15, 5)) -> None:
    """
    Compare kcat distributions between two datasets with comprehensive statistical analysis and visualization.
    
    Creates log-histogram with KDE overlay, ECDF plots, and Q-Q plot on log10 scale.
    Prints summary statistics for both original and log10-transformed values.
    
    Parameters
    ----------
    df1 : pd.DataFrame
        First dataframe containing kcat values
    kcat_col1 : str
        Column name containing kcat values in df1
    df2 : pd.DataFrame
        Second dataframe containing kcat values  
    kcat_col2 : str
        Column name containing kcat values in df2
    label1 : str, default "Dataset 1"
        Label for first dataset in plots and summary
    label2 : str, default "Dataset 2"
        Label for second dataset in plots and summary
    figsize : Tuple[int, int], default (15, 5)
        Figure size for the plot grid
    
    Returns
    -------
    None
        Function prints summary statistics and displays plots
    """
    
    # Extract and clean kcat values
    kcat1_raw = df1[kcat_col1].dropna()
    kcat2_raw = df2[kcat_col2].dropna()
    
    # Filter out non-positive values (required for log transformation)
    kcat1_raw = kcat1_raw[kcat1_raw > 0]
    kcat2_raw = kcat2_raw[kcat2_raw > 0]
    
    if len(kcat1_raw) == 0 or len(kcat2_raw) == 0:
        print("Error: No valid positive kcat values found in one or both datasets")
        return
    
    # Log10 transformation
    log_kcat1 = np.log10(kcat1_raw)
    log_kcat2 = np.log10(kcat2_raw)
    
    # Print summary statistics
    _print_summary_statistics(kcat1_raw, log_kcat1, kcat2_raw, log_kcat2, label1, label2)
    
    # Create visualization
    fig, axes = plt.subplots(1, 2, figsize=figsize)
    
    # Plot 1: Log-histogram + KDE overlay
    _plot_histogram_kde(log_kcat1, log_kcat2, axes[0], label1, label2)
    
    # Plot 2: ECDF overlay
    #_plot_ecdf(log_kcat1, log_kcat2, axes[1], label1, label2)
    
    # Plot 3: Q-Q plot
    _plot_qq(log_kcat1, log_kcat2, axes[1], label1, label2)
    
    plt.tight_layout()
    plt.show()


def _calculate_statistics(values: np.ndarray) -> Dict[str, float]:
    """Calculate comprehensive statistics for a dataset."""
    return {
        'count': len(values),
        'mean': np.mean(values),
        'median': np.median(values),
        'std': np.std(values, ddof=1),
        'p10': np.percentile(values, 10),
        'p25': np.percentile(values, 25),
        'p75': np.percentile(values, 75),
        'p90': np.percentile(values, 90),
        'iqr': np.percentile(values, 75) - np.percentile(values, 25)
    }


def _calculate_geometric_stats(values: np.ndarray) -> Dict[str, float]:
    """Calculate geometric mean and geometric standard deviation."""
    # Geometric mean = exp(mean(log(values))) = 10^(mean(log10(values)))
    log_values = np.log10(values)
    geom_mean = 10 ** np.mean(log_values)
    
    # Geometric SD = exp(std(log(values))) = 10^(std(log10(values)))
    geom_std = 10 ** np.std(log_values, ddof=1)
    
    return {
        'geometric_mean': geom_mean,
        'geometric_std': geom_std
    }


def _print_summary_statistics(kcat1_raw: pd.Series, log_kcat1: pd.Series,
                             kcat2_raw: pd.Series, log_kcat2: pd.Series,
                             label1: str, label2: str) -> None:
    """Print comprehensive summary statistics for both datasets."""
    
    print("="*80)
    print("KCAT DISTRIBUTION COMPARISON SUMMARY")
    print("="*80)
    
    # Calculate statistics
    stats1_raw = _calculate_statistics(kcat1_raw.values)
    stats1_log = _calculate_statistics(log_kcat1.values)
    geom_stats1 = _calculate_geometric_stats(kcat1_raw.values)
    
    stats2_raw = _calculate_statistics(kcat2_raw.values)
    stats2_log = _calculate_statistics(log_kcat2.values)
    geom_stats2 = _calculate_geometric_stats(kcat2_raw.values)
    
    # Print original scale statistics
    print("\nORIGINAL SCALE STATISTICS (s⁻¹)")
    print("-" * 50)
    print(f"{'Statistic':<20} {label1:<15} {label2:<15}")
    print("-" * 50)
    print(f"{'Count':<20} {stats1_raw['count']:<15,} {stats2_raw['count']:<15,}")
    print(f"{'Mean':<20} {stats1_raw['mean']:<15.2e} {stats2_raw['mean']:<15.2e}")
    print(f"{'Median':<20} {stats1_raw['median']:<15.2e} {stats2_raw['median']:<15.2e}")
    print(f"{'Std Dev':<20} {stats1_raw['std']:<15.2e} {stats2_raw['std']:<15.2e}")
    print(f"{'P10':<20} {stats1_raw['p10']:<15.2e} {stats2_raw['p10']:<15.2e}")
    print(f"{'P25':<20} {stats1_raw['p25']:<15.2e} {stats2_raw['p25']:<15.2e}")
    print(f"{'P75':<20} {stats1_raw['p75']:<15.2e} {stats2_raw['p75']:<15.2e}")
    print(f"{'P90':<20} {stats1_raw['p90']:<15.2e} {stats2_raw['p90']:<15.2e}")
    print(f"{'IQR (P25-P75)':<20} {stats1_raw['iqr']:<15.2e} {stats2_raw['iqr']:<15.2e}")
    print(f"{'Geometric Mean':<20} {geom_stats1['geometric_mean']:<15.2e} {geom_stats2['geometric_mean']:<15.2e}")
    print(f"{'Geometric Std':<20} {geom_stats1['geometric_std']:<15.2e} {geom_stats2['geometric_std']:<15.2e}")
    
    # Print log10 scale statistics
    print("\nLOG₁₀ SCALE STATISTICS")
    print("-" * 50)
    print(f"{'Statistic':<20} {label1:<15} {label2:<15}")
    print("-" * 50)
    print(f"{'Count':<20} {stats1_log['count']:<15,} {stats2_log['count']:<15,}")
    print(f"{'Mean':<20} {stats1_log['mean']:<15.3f} {stats2_log['mean']:<15.3f}")
    print(f"{'Median':<20} {stats1_log['median']:<15.3f} {stats2_log['median']:<15.3f}")
    print(f"{'Std Dev':<20} {stats1_log['std']:<15.3f} {stats2_log['std']:<15.3f}")
    print(f"{'P10':<20} {stats1_log['p10']:<15.3f} {stats2_log['p10']:<15.3f}")
    print(f"{'P25':<20} {stats1_log['p25']:<15.3f} {stats2_log['p25']:<15.3f}")
    print(f"{'P75':<20} {stats1_log['p75']:<15.3f} {stats2_log['p75']:<15.3f}")
    print(f"{'P90':<20} {stats1_log['p90']:<15.3f} {stats2_log['p90']:<15.3f}")
    print(f"{'IQR (P25-P75)':<20} {stats1_log['iqr']:<15.3f} {stats2_log['iqr']:<15.3f}")
    
    print("="*80)


def _plot_histogram_kde(log_kcat1: pd.Series, log_kcat2: pd.Series, ax: plt.Axes,
                       label1: str, label2: str) -> None:
    """Plot log-histogram with KDE overlay."""
    
    # Determine appropriate number of bins
    n_bins = max(20, min(50, int(np.sqrt(len(log_kcat1) + len(log_kcat2)))))
    
    # Plot histograms
    ax.hist(log_kcat1, bins=n_bins, alpha=0.6, density=True, 
            label=f'{label1} (n={len(log_kcat1):,})', color='skyblue', edgecolor='black', linewidth=0.5)
    ax.hist(log_kcat2, bins=n_bins, alpha=0.6, density=True,
            label=f'{label2} (n={len(log_kcat2):,})', color='lightcoral', edgecolor='black', linewidth=0.5)
    
    # Add KDE overlays
    try:
        kde1 = stats.gaussian_kde(log_kcat1.dropna())
        kde2 = stats.gaussian_kde(log_kcat2.dropna())
        
        x_range = np.linspace(min(log_kcat1.min(), log_kcat2.min()),
                             max(log_kcat1.max(), log_kcat2.max()), 200)
        
        ax.plot(x_range, kde1(x_range), color='blue', linewidth=2, linestyle='--')
        ax.plot(x_range, kde2(x_range), color='red', linewidth=2, linestyle='--')
        
    except Exception as e:
        print(f"Warning: Could not generate KDE overlay: {e}")
    
    ax.set_xlabel('log₁₀(kcat) [s⁻¹]', fontweight="bold")
    ax.set_ylabel('Density', fontweight="bold")
    ax.set_title('Distribution Comparison\n(Histogram + KDE)')
    ax.legend()
    ax.grid(True, alpha=0.3)


def _plot_ecdf(log_kcat1: pd.Series, log_kcat2: pd.Series, ax: plt.Axes,
               label1: str, label2: str) -> None:
    """Plot empirical cumulative distribution functions."""
    
    def ecdf(data):
        """Calculate empirical CDF."""
        x = np.sort(data)
        y = np.arange(1, len(x) + 1) / len(x)
        return x, y
    
    x1, y1 = ecdf(log_kcat1.dropna())
    x2, y2 = ecdf(log_kcat2.dropna())
    
    ax.plot(x1, y1, label=f'{label1} (n={len(log_kcat1):,})', 
            linewidth=2, color='blue')
    ax.plot(x2, y2, label=f'{label2} (n={len(log_kcat2):,})', 
            linewidth=2, color='red')
    
    ax.set_xlabel('log₁₀(kcat) [s⁻¹]', fontweight="bold")
    ax.set_ylabel('Cumulative Probability', fontweight="bold")
    ax.set_title('Empirical Cumulative\nDistribution Functions')
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_ylim(0, 1)



def _plot_qq(log_kcat_x, log_kcat_y, ax: plt.Axes,
             label_x: str, label_y: str,
             qmin: float = 0.02, qmax: float = 0.98, n_q: int = 49,
             random_state: int = 0) -> None:
    """
    Quantile–Quantile plot on log10(kcat):
      x-axis: label_x (e.g., kcat_invivo)
      y-axis: label_y (e.g., kcat_CatPred)

    Shows y=x (no distributional difference) and a robust Theil–Sen fit:
        q_y(p) ≈ a + b * q_x(p)
    where:
      a  ~ location shift on log10 scale (10**a is the multiplicative shift)
      b  ~ relative dispersion (b≈1 similar spread; b<1 narrower; b>1 wider)
    """

    # clean & choose quantiles
    x = np.asarray(log_kcat_x.dropna(), dtype=float)
    y = np.asarray(log_kcat_y.dropna(), dtype=float)
    if x.size == 0 or y.size == 0:
        raise ValueError("Empty input after dropping NaNs.")

    p = np.linspace(qmin, qmax, min(n_q, x.size, y.size))
    qx = np.quantile(x, p)
    qy = np.quantile(y, p)

    # scatter of Q–Q points
    ax.scatter(qx, qy, alpha=0.65, s=22, label="Quantile pairs", color="purple")

    # y = x reference line
    lo = float(min(qx.min(), qy.min()))
    hi = float(max(qx.max(), qy.max()))
    ax.plot([lo, hi], [lo, hi], "r--", lw=2, label="y = x")

    # ordinary least squares fit
    a = b = np.nan
    b, a = np.polyfit(qx, qy, deg=1)

    # plot the fitted line across the Q–Q domain
    xs = np.array([lo, hi])
    ax.plot(xs, a + b * xs, lw=2.0, color="black",
            label=("y = a + b x"))

    # annotate interpretable stats
    mad_resid = np.median(np.abs(qy - (a + b * qx)))  # robust residual scale on log10
    ax.text(0.03, 0.97,
            f"a = {a:.2f}\n"
            f"b = {b:.2f}\n"
            f"MAD = {mad_resid:.2f}",
            transform=ax.transAxes, ha="left", va="top",
            fontsize=10, bbox=dict(boxstyle="round", facecolor="white", alpha=0.85))

    # axes/labels 
    ax.set_xlabel("log₁₀(kcat) quantiles — in vivo", fontweight="bold")
    ax.set_ylabel("log₁₀(kcat) quantiles — in vitro", fontweight="bold")
    ax.set_title("Quantile–Quantile Plot (log₁₀ scale)")
    ax.grid(True, alpha=0.3)
    ax.legend(frameon=True)
    



# ============================================
# Functions for specific kcat datasets
# ============================================

def load_kcat_dataset_ecoli(CPIPred_dir, CatPred_dir, EnzyExtract_dir) -> pd.DataFrame:
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



# ============================================
# Others
# ============================================

def canonicalize(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        return Chem.MolToSmiles(mol, canonical=True)
    except:
        return None