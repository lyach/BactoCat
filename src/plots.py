"""
Module Description: plots.py

Purpose: 
Utility functions for kcat visualization.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
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


def _plot_histogram_kde(log_kcat1: pd.Series, log_kcat2: pd.Series = None, ax: plt.Axes = None,
                       label1: str = "Dataset", label2: str = "Dataset 2") -> None:
    """
    Plot log-histogram with KDE overlay for one or two datasets.
    
    Parameters
    ----------
    log_kcat1 : pd.Series
        First dataset (log-transformed kcat values)
    log_kcat2 : pd.Series, optional
        Second dataset for comparison. If None, only plots log_kcat1
    ax : plt.Axes, optional
        Matplotlib axes to plot on. If None, creates new figure
    label1 : str, default "Dataset"
        Label for first dataset
    label2 : str, default "Dataset 2"
        Label for second dataset (only used if log_kcat2 is provided)
    """
    
    # Create figure if axes not provided
    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 6))
    
    # Determine if plotting one or two datasets
    dual_mode = log_kcat2 is not None
    
    # Number of bins
    if dual_mode:
        n_bins = max(20, min(50, int(np.sqrt(len(log_kcat1) + len(log_kcat2)))))
    else:
        n_bins = max(20, min(50, int(np.sqrt(len(log_kcat1)))))
    
    # Plot histogram(s)
    if dual_mode:
        ax.hist(log_kcat1, bins=n_bins, alpha=0.6, density=True, 
                label=f'{label1} (n={len(log_kcat1):,})', color='skyblue', edgecolor='black', linewidth=0.5)
        ax.hist(log_kcat2, bins=n_bins, alpha=0.6, density=True,
                label=f'{label2} (n={len(log_kcat2):,})', color='lightcoral', edgecolor='black', linewidth=0.5)
    else:
        ax.hist(log_kcat1, bins=n_bins, alpha=0.7, density=True, 
                label=f'{label1} (n={len(log_kcat1):,})', color='skyblue', edgecolor='black', linewidth=0.5)
    
    # Add KDE overlay
    try:
        kde1 = stats.gaussian_kde(log_kcat1.dropna())
        x_min = log_kcat1.min()
        x_max = log_kcat1.max()
        
        if dual_mode:
            kde2 = stats.gaussian_kde(log_kcat2.dropna())
            x_min = min(x_min, log_kcat2.min())
            x_max = max(x_max, log_kcat2.max())
            
            x_range = np.linspace(x_min, x_max, 200)
            ax.plot(x_range, kde1(x_range), color='blue', linewidth=2, linestyle='--')
            ax.plot(x_range, kde2(x_range), color='red', linewidth=2, linestyle='--')
        else:
            x_range = np.linspace(x_min, x_max, 200)
            ax.plot(x_range, kde1(x_range), color='blue', linewidth=2, linestyle='--')
        
    except Exception as e:
        print(f"Warning: Could not generate KDE overlay: {e}")
    
    ax.set_xlabel('log₁₀(kcat) [s⁻¹]', fontweight="bold")
    ax.set_ylabel('Density', fontweight="bold")
    
    if dual_mode:
        ax.set_title('Distribution Comparison\n(Histogram + KDE)')
    else:
        ax.set_title('kcat Distribution\n(Histogram + KDE)')
    
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    if ax is None:
        plt.tight_layout()
        plt.show()


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
    

def plot_eta_variability(df: pd.DataFrame, figsize: Tuple[int, int] = (12, 3)):
    """
    Visualize eta variability across enzyme-substrate pairs with comprehensive plots.
    
    Creates a 3 panel figure showing:
    - Distribution of eta_mean values
    - Coefficient of variation (CV) analysis
    - Scatter plots of variance metrics
    - Range analysis (min vs max)
    
    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing eta variance metrics with columns:
        'eta_mean', 'eta_stdev', 'eta_min', 'eta_max', 'eta_cv'
    figsize : Tuple[int, int], default (18, 12)
        Figure size for the plot grid
    
    Returns
    -------
    None
        Function displays plots and prints summary statistics
    """
    
    # Check required columns
    required_cols = ['eta_mean', 'eta_stdev', 'eta_min', 'eta_max', 'eta_cv']
    missing_cols = [col for col in required_cols if col not in df.columns]
    
    # Filter out rows with NaN values in eta_mean
    df_clean = df[df['eta_mean'].notna()].copy()
    
    
    print("="*80)
    print("ETA VARIABILITY ANALYSIS")
    print("="*80)
    print(f"Total enzyme-substrate pairs analyzed: {len(df_clean):,}")
    
    print(f"\nEta Coefficient of Variation (CV) Statistics:")
    cv_clean = df_clean[df_clean['eta_cv'].notna()]
    if len(cv_clean) > 0:
        print(f"  Mean:   {cv_clean['eta_cv'].mean():.3f}")
        print(f"  Median: {cv_clean['eta_cv'].median():.3f}")
        print(f"  Min:    {cv_clean['eta_cv'].min():.3f}")
        print(f"  Max:    {cv_clean['eta_cv'].max():.3f}")
    else:
        print("  No valid CV values")
    
    print(f"\nEta Range Statistics:")
    print(f"  Min eta across all conditions: {df_clean['eta_min'].min():.3f}")
    print(f"  Max eta across all conditions: {df_clean['eta_max'].max():.3f}")
    print("="*80)
    
    # Create the visualization
    fig = plt.figure(figsize=figsize)
    gs = fig.add_gridspec(1, 2, hspace=0.3, wspace=0.3)
    
    # Plot 1: Distribution of eta_mean 
    ax1 = fig.add_subplot(gs[0, 0])
    _plot_eta_mean_distribution(df_clean, ax1)
    
    # Plot 2: Distribution of eta_cv 
    ax2 = fig.add_subplot(gs[0, 1])
    _plot_eta_cv_distribution(df_clean, ax2)
    
    # Plot 3: Eta mean vs CV scatter
    # ax3 = fig.add_subplot(gs[0, 2])
    # _plot_eta_mean_vs_cv(df_clean, ax3)
    
    plt.show()


def _plot_eta_mean_distribution(df: pd.DataFrame, ax: plt.Axes):
    """Plot distribution of eta_mean values."""
    data = df['eta_mean'].dropna()
    
    # Histogram
    ax.hist(data, bins=30, alpha=0.7, color='steelblue', edgecolor='black', 
            density=True, label=f'n={len(data):,}')
    
    # KDE overlay
    try:
        kde = stats.gaussian_kde(data)
        x_range = np.linspace(data.min(), data.max(), 200)
        ax.plot(x_range, kde(x_range), 'r-', linewidth=2)
    except:
        pass
    
    # Add vertical line at mean
    ax.axvline(data.mean(), color='red', linestyle='--', linewidth=2, 
               label=f'Mean={data.mean():.3f}')
    ax.axvline(data.median(), color='green', linestyle='--', linewidth=2,
               label=f'Median={data.median():.3f}')
    
    ax.set_xlabel('η mean', fontweight='bold')
    ax.set_ylabel('Density', fontweight='bold')
    ax.set_title('Distribution of Mean η Values')
    ax.legend()
    ax.grid(True, alpha=0.3)


def _plot_eta_cv_distribution(df: pd.DataFrame, ax: plt.Axes):
    """Plot distribution of coefficient of variation."""
    data = df['eta_cv'].dropna()
    
    if len(data) == 0:
        ax.text(0.5, 0.5, 'No CV data available', ha='center', va='center',
                transform=ax.transAxes)
        return
    
    # Histogram
    ax.hist(data, bins=30, alpha=0.7, color='coral', edgecolor='black',
            density=True, label=f'n={len(data):,}')
    
    # KDE overlay
    try:
        kde = stats.gaussian_kde(data)
        x_range = np.linspace(data.min(), data.max(), 200)
        ax.plot(x_range, kde(x_range), 'darkred', linewidth=2)
    except:
        pass
    
    ax.axvline(data.mean(), color='red', linestyle='--', linewidth=2,
               label=f'Mean={data.mean():.3f}')
    
    ax.set_xlabel('Coefficient of Variation (CV)', fontweight='bold')
    ax.set_ylabel('Density', fontweight='bold')
    ax.set_title('Distribution of η CV')
    ax.legend()
    ax.grid(True, alpha=0.3)


def _plot_eta_mean_vs_cv(df: pd.DataFrame, ax: plt.Axes):
    """Scatter plot of eta_mean vs eta_cv."""
    data = df[['eta_mean', 'eta_cv']].dropna()
    
    if len(data) == 0:
        ax.text(0.5, 0.5, 'No data available', ha='center', va='center',
                transform=ax.transAxes)
        return
    
    # Scatter plot with alpha for overlapping points
    scatter = ax.scatter(data['eta_mean'], data['eta_cv'], 
                        alpha=0.5, s=20, c=data['eta_cv'], 
                        cmap='viridis', edgecolors='k', linewidth=0.5)
    
    # Add colorbar
    plt.colorbar(scatter, ax=ax, label='CV')
    
    # Add trend line if enough points
    if len(data) > 10:
        z = np.polyfit(data['eta_mean'], data['eta_cv'], 1)
        p = np.poly1d(z)
        x_trend = np.linspace(data['eta_mean'].min(), data['eta_mean'].max(), 100)
        ax.plot(x_trend, p(x_trend), 'r--', linewidth=2, alpha=0.7,
                label=f'y={z[0]:.2f}x+{z[1]:.2f}')
    
    ax.set_xlabel('η mean', fontweight='bold')
    ax.set_ylabel('Coefficient of Variation (CV)', fontweight='bold')
    ax.set_title('η Mean vs Variability')
    ax.grid(True, alpha=0.3)
    if len(data) > 10:
        ax.legend()

def group_eta_variability(csv_path, deduplicate=False):
    """
    Load in vivo kcat with eta and stratify into high vs low variance groups.
    
    Parameters
    ----------
    csv_path : str or Path
        Path to the CSV file containing the in vivo variability data
        
    Returns
    -------
    pd.DataFrame
        Deduplicated dataframe with variance_group column added
    """
    # Load dataset
    df = pd.read_csv(csv_path)
    
    # Stratify into high vs low variance groups
    high_var_threshold = df['eta_cv'].quantile(0.75)
    low_var_threshold = df['eta_cv'].quantile(0.25)
    
    # Groups: medium, high, low
    df['variance_group'] = 'medium'
    df.loc[df['eta_cv'] >= high_var_threshold, 'variance_group'] = 'high'
    df.loc[df['eta_cv'] <= low_var_threshold, 'variance_group'] = 'low'
    
    if deduplicate:
        # Deduplicate by gene (keeping first)
        print(f"Total rows before deduplication: {len(df)}")
        df_unique = df.drop_duplicates(subset='gene', keep='first')
        print(f"Total rows after deduplication: {len(df_unique)}")
        print(f"{len(df) - len(df_unique)} duplicate genes removed")
        df = df_unique
    
    # Show breakdown by group
    print(f"\nHigh variance: {(df['variance_group']=='high').sum()} enzymes")
    print(f"Low variance: {(df['variance_group']=='low').sum()} enzymes")
    print(f"Medium variance: {(df['variance_group']=='medium').sum()} enzymes")
    
    return df
