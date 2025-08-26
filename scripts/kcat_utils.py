"""
Module Description: kcat_utils.py

Purpose: 
Utility functions for kcat visualization.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from typing import Tuple, Dict, Any
import warnings
warnings.filterwarnings('ignore', category=RuntimeWarning)


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
    fig, axes = plt.subplots(1, 3, figsize=figsize)
    
    # Plot 1: Log-histogram + KDE overlay
    _plot_histogram_kde(log_kcat1, log_kcat2, axes[0], label1, label2)
    
    # Plot 2: ECDF overlay
    _plot_ecdf(log_kcat1, log_kcat2, axes[1], label1, label2)
    
    # Plot 3: Q-Q plot
    _plot_qq(log_kcat1, log_kcat2, axes[2], label1, label2)
    
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
    print(f"\nORIGINAL SCALE STATISTICS (s⁻¹)")
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
    print(f"\nLOG₁₀ SCALE STATISTICS")
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
    
    ax.set_xlabel('log₁₀(kcat) [s⁻¹]')
    ax.set_ylabel('Density')
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
    
    ax.set_xlabel('log₁₀(kcat) [s⁻¹]')
    ax.set_ylabel('Cumulative Probability')
    ax.set_title('Empirical Cumulative\nDistribution Functions')
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_ylim(0, 1)


def _plot_qq(log_kcat1: pd.Series, log_kcat2: pd.Series, ax: plt.Axes,
             label1: str, label2: str) -> None:
    """Plot quantile-quantile plot."""
    
    # Clean data
    data1 = log_kcat1.dropna()
    data2 = log_kcat2.dropna()
    
    # Calculate quantiles for Q-Q plot
    n_quantiles = min(len(data1), len(data2), 1000)  # Limit for performance
    quantiles = np.linspace(0, 1, n_quantiles)
    
    q1 = np.quantile(data1, quantiles)
    q2 = np.quantile(data2, quantiles)
    
    # Plot Q-Q
    ax.scatter(q1, q2, alpha=0.6, s=20, color='purple')
    
    # Add diagonal reference line
    min_val = min(q1.min(), q2.min())
    max_val = max(q1.max(), q2.max())
    ax.plot([min_val, max_val], [min_val, max_val], 
            'r--', linewidth=2, label='y = x')
    
    ax.set_xlabel(f'log₁₀(kcat) Quantiles - {label1}')
    ax.set_ylabel(f'log₁₀(kcat) Quantiles - {label2}')
    ax.set_title('Quantile-Quantile Plot\n(log₁₀ scale)')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Add correlation coefficient
    correlation = np.corrcoef(q1, q2)[0, 1]
    ax.text(0.05, 0.95, f'r = {correlation:.3f}', 
            transform=ax.transAxes, fontsize=10,
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

