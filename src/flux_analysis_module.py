import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scipy.stats import pearsonr, spearmanr
import warnings
warnings.filterwarnings('ignore')

class FluxAnalyzer:
    """
    A comprehensive module for analyzing and comparing metabolic flux values
    from different methods (experimental, reconciliation, FBA, etc.)
    """
    
    def __init__(self, df, method1_col=None, method2_col=None, rxn_id_col=None):
        """
        Initialize the FluxAnalyzer with a dataframe
        
        Parameters:
        -----------
        df : pandas.DataFrame
            DataFrame with 3 columns: reaction_id, flux_method1, flux_method2
        method1_col : str, optional
            Name of first flux column (auto-detected if None)
        method2_col : str, optional
            Name of second flux column (auto-detected if None)
        rxn_id_col : str, optional
            Name of reaction ID column (auto-detected if None)
        """
        self.df = df.copy()
        
        # Auto-detect columns if not specified
        if rxn_id_col is None:
            rxn_id_col = df.columns[0]
        if method1_col is None:
            method1_col = df.columns[1]
        if method2_col is None:
            method2_col = df.columns[2]
            
        self.rxn_id_col = rxn_id_col
        self.method1_col = method1_col
        self.method2_col = method2_col
        
        # Clean and prepare data
        self._prepare_data()
        
    def _prepare_data(self):
        """Clean and prepare data for analysis"""
        # Convert flux columns to numeric, handling any non-numeric values
        self.df[self.method1_col] = pd.to_numeric(self.df[self.method1_col], errors='coerce')
        self.df[self.method2_col] = pd.to_numeric(self.df[self.method2_col], errors='coerce')
        
        # Create a clean dataset with no NaN values for correlation analysis
        self.clean_df = self.df.dropna(subset=[self.method1_col, self.method2_col])
        
        # Calculate absolute differences
        self.clean_df = self.clean_df.copy()
        self.clean_df['abs_difference'] = np.abs(self.clean_df[self.method1_col] - self.clean_df[self.method2_col])
        self.clean_df['relative_difference'] = self.clean_df['abs_difference'] / (np.abs(self.clean_df[self.method1_col]) + 1e-10)
        
    def correlation_analysis(self, plot=True, figsize=(10, 6)):
        """
        Perform correlation analysis between two flux methods
        
        Parameters:
        -----------
        plot : bool
            Whether to create correlation plots
        figsize : tuple
            Figure size for plots
            
        Returns:
        --------
        dict : Dictionary containing correlation statistics
        """
        if len(self.clean_df) == 0:
            print("No valid data points for correlation analysis")
            return {}
            
        # Calculate correlations
        pearson_r, pearson_p = pearsonr(self.clean_df[self.method1_col], self.clean_df[self.method2_col])
        spearman_r, spearman_p = spearmanr(self.clean_df[self.method1_col], self.clean_df[self.method2_col])
        
        # R-squared
        r2 = pearson_r ** 2
        
        # Linear regression
        slope, intercept, _, _, _ = stats.linregress(self.clean_df[self.method1_col], self.clean_df[self.method2_col])
        
        results = {
            'pearson_r': pearson_r,
            'pearson_p': pearson_p,
            'spearman_r': spearman_r,
            'spearman_p': spearman_p,
            'r_squared': r2,
            'slope': slope,
            'intercept': intercept,
            'n_points': len(self.clean_df)
        }
        
        if plot:
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)
            
            # Scatter plot with regression line
            ax1.scatter(self.clean_df[self.method1_col], self.clean_df[self.method2_col], 
                       alpha=0.7, s=60, edgecolors='black', linewidth=0.5)
            
            # Add regression line
            x_range = np.linspace(self.clean_df[self.method1_col].min(), 
                                 self.clean_df[self.method1_col].max(), 100)
            y_pred = slope * x_range + intercept
            ax1.plot(x_range, y_pred, 'r--', linewidth=2, label=f'y = {slope:.3f}x + {intercept:.3f}')
            
            # Add diagonal line (perfect correlation)
            min_val = min(self.clean_df[self.method1_col].min(), self.clean_df[self.method2_col].min())
            max_val = max(self.clean_df[self.method1_col].max(), self.clean_df[self.method2_col].max())
            ax1.plot([min_val, max_val], [min_val, max_val], 'k--', alpha=0.5, label='Perfect correlation')
            
            ax1.set_xlabel(f'{self.method1_col}')
            ax1.set_ylabel(f'{self.method2_col}')
            ax1.set_title(f'Correlation Plot\nR² = {r2:.3f}, r = {pearson_r:.3f}')
            ax1.legend()
            ax1.grid(True, alpha=0.3)
            
            # Residual plot
            residuals = self.clean_df[self.method2_col] - (slope * self.clean_df[self.method1_col] + intercept)
            ax2.scatter(self.clean_df[self.method1_col], residuals, alpha=0.7, s=60, 
                       edgecolors='black', linewidth=0.5)
            ax2.axhline(y=0, color='r', linestyle='--', linewidth=2)
            ax2.set_xlabel(f'{self.method1_col}')
            ax2.set_ylabel('Residuals')
            ax2.set_title('Residual Plot')
            ax2.grid(True, alpha=0.3)
            
            plt.tight_layout()
            plt.show()
            
        return results
    
    def plot_correlation_enhanced(self, figsize=(8, 6), sizes=(20, 200), palette='coolwarm', alpha=0.7, 
                                 xlabel=None, ylabel=None):
        """
        Create an enhanced correlation plot with dots colored and scaled by absolute difference
        
        Parameters:
        -----------
        figsize : tuple
            Figure size for the plot
        sizes : tuple
            Min and max sizes for the scatter points
        palette : str
            Color palette for the hue mapping
        alpha : float
            Transparency of the scatter points
        xlabel : str, optional
            Custom label for x-axis. If None, uses column name
        ylabel : str, optional
            Custom label for y-axis. If None, uses column name
            
        Returns:
        --------
        dict : Dictionary containing correlation statistics
        """
        if len(self.clean_df) == 0:
            print("No valid data points for correlation analysis")
            return {}
        
        # Calculate correlation statistics
        pearson_r, pearson_p = pearsonr(self.clean_df[self.method1_col], self.clean_df[self.method2_col])
        r2 = pearson_r ** 2
        
        # Linear regression
        slope, intercept, _, _, _ = stats.linregress(self.clean_df[self.method1_col], self.clean_df[self.method2_col])
        
        # Create the plot
        plt.figure(figsize=figsize)
        
        # Create enhanced scatter plot with size and color based on absolute difference
        sns.scatterplot(data=self.clean_df, x=self.method1_col, y=self.method2_col,
                       size='abs_difference', hue='abs_difference',
                       sizes=sizes, palette=palette, alpha=alpha)
        
        # Add regression line
        x_range = np.linspace(self.clean_df[self.method1_col].min(), 
                             self.clean_df[self.method1_col].max(), 100)
        y_pred = slope * x_range + intercept
        plt.plot(x_range, y_pred, 'r--', linewidth=2, label=f'y = {slope:.3f}x + {intercept:.3f}')
        
        # Fix the legend title
        legend = plt.legend(title='Abs difference', loc='upper left', 
                           bbox_to_anchor=(0.02, 0.98),
                           title_fontsize=11, fontsize=10,
                           framealpha=0.9, fancybox=True)
        
        # Add perfect correlation line
        min_val = min(self.clean_df[self.method1_col].min(), self.clean_df[self.method2_col].min())
        max_val = max(self.clean_df[self.method1_col].max(), self.clean_df[self.method2_col].max())
        plt.plot([min_val, max_val], [min_val, max_val], 'k--', alpha=0.7, label='Perfect correlation')
        
        # Set labels and title
        plt.xlabel(xlabel if xlabel is not None else f'{self.method1_col}')
        plt.ylabel(ylabel if ylabel is not None else f'{self.method2_col}')
        plt.title(f'Enhanced Correlation Plot\nR² = {r2:.3f}, r = {pearson_r:.3f}')
        plt.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.show()
        
        # Return correlation statistics
        results = {
            'pearson_r': pearson_r,
            'pearson_p': pearson_p,
            'r_squared': r2,
            'slope': slope,
            'intercept': intercept,
            'n_points': len(self.clean_df)
        }
        
        return results
    
    def identify_outliers(self, top_n=10, method='absolute'):
        """
        Identify reactions with biggest and smallest differences
        
        Parameters:
        -----------
        top_n : int
            Number of top reactions to return
        method : str
            'absolute' for absolute differences, 'relative' for relative differences
            
        Returns:
        --------
        dict : Dictionary with biggest and smallest differences, plus the full sorted dataframe
            - 'biggest_differences': Top N reactions with largest differences
            - 'smallest_differences': Top N reactions with smallest differences  
            - 'sorted_df': Full dataframe sorted by difference (descending)
        """
        if method == 'absolute':
            diff_col = 'abs_difference'
        elif method == 'relative':
            diff_col = 'relative_difference'
        else:
            raise ValueError("Method must be 'absolute' or 'relative'")
            
        # Sort by difference
        sorted_df = self.clean_df.sort_values(diff_col, ascending=False)
        
        biggest_diff = sorted_df.head(top_n)[[self.rxn_id_col, self.method1_col, self.method2_col, diff_col]]
        smallest_diff = sorted_df.tail(top_n)[[self.rxn_id_col, self.method1_col, self.method2_col, diff_col]]
        
        return {
            'biggest_differences': biggest_diff,
            'smallest_differences': smallest_diff,
            'sorted_df': sorted_df
        }
    
    def statistical_summary(self):
        """
        Generate statistical summary of the flux comparisons
        
        Returns:
        --------
        dict : Dictionary with statistical measures
        """
        summary = {
            'total_reactions': len(self.df),
            'reactions_with_both_values': len(self.clean_df),
            'reactions_missing_method1': self.df[self.method1_col].isna().sum(),
            'reactions_missing_method2': self.df[self.method2_col].isna().sum(),
            'mean_absolute_difference': self.clean_df['abs_difference'].mean(),
            'median_absolute_difference': self.clean_df['abs_difference'].median(),
            'std_absolute_difference': self.clean_df['abs_difference'].std(),
            'mean_relative_difference': self.clean_df['relative_difference'].mean(),
            'median_relative_difference': self.clean_df['relative_difference'].median(),
            'max_absolute_difference': self.clean_df['abs_difference'].max(),
            'min_absolute_difference': self.clean_df['abs_difference'].min(),
        }
        
        return summary
    
    
    def generate_report(self, output_file=None):
        """
        Generate a comprehensive analysis report
        
        Parameters:
        -----------
        output_file : str, optional
            Path to save the report as text file
        """
        report = []
        report.append("="*60)
        report.append("METABOLIC FLUX ANALYSIS REPORT")
        report.append("="*60)
        report.append(f"Comparing: {self.method1_col} vs {self.method2_col}")
        report.append("")
        
        # Statistical summary
        summary = self.statistical_summary()
        report.append("STATISTICAL SUMMARY:")
        report.append("-" * 20)
        for key, value in summary.items():
            report.append(f"{key}: {value:.6f}" if isinstance(value, float) else f"{key}: {value}")
        report.append("")
        
        # Correlation analysis
        corr_results = self.correlation_analysis(plot=False)
        report.append("CORRELATION ANALYSIS:")
        report.append("-" * 20)
        for key, value in corr_results.items():
            report.append(f"{key}: {value:.6f}" if isinstance(value, float) else f"{key}: {value}")
        report.append("")
        
        # Outliers
        outliers = self.identify_outliers(top_n=10, method='absolute')
        report.append(f"TOP REACTIONS WITH BIGGEST ABSOLUTE DIFFERENCES ({self.method1_col} vs {self.method2_col}):")
        report.append("-" * 50)
        for _, row in outliers['biggest_differences'].iterrows():
            report.append(f"{row[self.rxn_id_col]}: {row[self.method1_col]:.6f} vs {row[self.method2_col]:.6f} (diff: {row['abs_difference']:.6f})")
        report.append("")
        
        report.append(f"TOP REACTIONS WITH SMALLEST ABSOLUTE DIFFERENCES ({self.method1_col} vs {self.method2_col}):")
        report.append("-" * 50)
        for _, row in outliers['smallest_differences'].iterrows():
            report.append(f"{row[self.rxn_id_col]}: {row[self.method1_col]:.6f} vs {row[self.method2_col]:.6f} (diff: {row['abs_difference']:.6f})")
        
        report_text = "\n".join(report)
        
        if output_file:
            with open(output_file, 'w') as f:
                f.write(report_text)
            print(f"Report saved to {output_file}")
        else:
            print(report_text)
        
        return report_text


# Example of how to use with your data
if __name__ == "__main__":
    pass