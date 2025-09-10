
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scipy.stats import wilcoxon, ks_2samp
import cobra
from cobra.io import read_sbml_model
from cobra.flux_analysis import pfba


def load_model(model_id):
    """
    Load a metabolic model using COBRA.
    
    Parameters:
    -----------
    model_id : str
        ID of the model (e.g., 'iML1515') to load using COBRA's model repository,
        or path to the model file if not using the cobra repository.
        
    Returns:
    --------
    cobra.Model or None
        The loaded COBRA model, or None if loading failed.
    """
    try:
        # From the cobra repository
        model = cobra.io.load_model(model_id)
        print(f"Model {model_id} loaded successfully from COBRA repository.")
        return model
    except Exception:
        # Try alternate loading method (model ID is a file path)
        print(f"Loading model as path...")
        try:
            model = cobra.io.read_sbml_model(model_id)
            print(f"Model {model_id} loaded successfully from path.")
            return model
        except Exception as e2:
            print(f"Failed to load model: {e2}")
            return None
        
        
def compare_pFBA_FBA(model_id: str):
    '''
    Compare pFBA and FBA solutions.
    
    Parameters:
        model: cobra.Model
            The model to compare.
            
    Returns:
        df_fluxes: pd.DataFrame
    '''
    
    model = load_model(model_id)
    
    # ----- Run optimizations ------------------------------------------------
    solution_fba  = model.optimize()
    solution_pfba = pfba(model)

    # ----- Build dataframe --------------------------------------------------
    df_fluxes = pd.DataFrame({
        "reaction": solution_pfba.fluxes.index,
        "pFBA":    solution_pfba.fluxes.values,
        "FBA":     solution_fba.fluxes.values,
    })
    df_fluxes["abs_pFBA"]       = np.abs(df_fluxes["pFBA"])
    df_fluxes["abs_FBA"]        = np.abs(df_fluxes["FBA"])
    df_fluxes["difference"]     = df_fluxes["FBA"] - df_fluxes["pFBA"]
    df_fluxes["abs_difference"] = np.abs(df_fluxes["difference"])
    df_fluxes["active_pFBA"]    = df_fluxes["abs_pFBA"] > 1e-6
    df_fluxes["active_FBA"]     = df_fluxes["abs_FBA"]  > 1e-6

    # -----------------------------------------------------------------------
    # 1. BASIC STATS
    # -----------------------------------------------------------------------
    print("\n=== CRITICAL COMPARISON OF FBA vs pFBA ===")
    print("1. BASIC STATISTICS")
    print(f"Total reactions:                 {len(df_fluxes)}")
    print(f"Active reactions (FBA):          {df_fluxes['active_FBA'].sum()}")
    print(f"Active reactions (pFBA):         {df_fluxes['active_pFBA'].sum()}")
    print(f"Sum |flux| (FBA):                {df_fluxes['abs_FBA'].sum():.4f}")
    print(f"Sum |flux| (pFBA):               {df_fluxes['abs_pFBA'].sum():.4f}")
    print(f"Objective value (FBA):           {solution_fba.objective_value:.4f}")
    print(f"Objective value (pFBA):          {solution_pfba.objective_value:.4f}")

    # -----------------------------------------------------------------------
    # 2. CORRELATION
    # -----------------------------------------------------------------------
    print("\n2. CORRELATION ANALYSIS")
    pearson_r, pearson_p   = stats.pearsonr(df_fluxes["FBA"], df_fluxes["pFBA"])
    spearman_r, spearman_p = stats.spearmanr(df_fluxes["FBA"], df_fluxes["pFBA"])
    print(f"Pearson  r={pearson_r:.4f}, p={pearson_p:.2e}")
    print(f"Spearman ρ={spearman_r:.4f}, p={spearman_p:.2e}")

    # -----------------------------------------------------------------------
    # 3. FLUX MAGNITUDE TEST (Wilcoxon)
    # -----------------------------------------------------------------------
    print("\n3. FLUX MAGNITUDE COMPARISON")
    active_mask   = df_fluxes["active_FBA"] | df_fluxes["active_pFBA"]
    active_fluxes = df_fluxes[active_mask]
    if len(active_fluxes):
        wil_stat, wil_p = wilcoxon(
            active_fluxes["abs_FBA"], active_fluxes["abs_pFBA"], alternative="two-sided"
        )
        interp = "Significant" if wil_p < 0.05 else "Not significant"
        print(f"Wilcoxon W={wil_stat:.4f}, p={wil_p:.2e} → {interp}")

    # -----------------------------------------------------------------------
    # 4. ACTIVE-/INACTIVE REACTION CONTINGENCY (McNemar)
    # -----------------------------------------------------------------------
    print("\n4. ACTIVE REACTION COMPARISON")
    both_active  = (df_fluxes["active_FBA"] &  df_fluxes["active_pFBA"]).sum()
    fba_only     = (df_fluxes["active_FBA"] & ~df_fluxes["active_pFBA"]).sum()
    pfba_only    = (~df_fluxes["active_FBA"] &  df_fluxes["active_pFBA"]).sum()
    print(f"Both active   : {both_active}")
    print(f"FBA-only      : {fba_only}")
    print(f"pFBA-only     : {pfba_only}")

    if fba_only + pfba_only:
        m_stat = (abs(fba_only - pfba_only) - 1) ** 2 / (fba_only + pfba_only)
        m_p    = 1 - stats.chi2.cdf(m_stat, 1)
        print(f"McNemar χ²={m_stat:.4f}, p={m_p:.4f}")

    # -----------------------------------------------------------------------
    # 5. DISTRIBUTION TEST (Kolmogorov-Smirnov)
    # -----------------------------------------------------------------------
    print("\n5. DISTRIBUTION COMPARISON")
    nonzero_fba  = df_fluxes.loc[df_fluxes["active_FBA"],  "FBA"]
    nonzero_pfba = df_fluxes.loc[df_fluxes["active_pFBA"], "pFBA"]
    if len(nonzero_fba) and len(nonzero_pfba):
        ks_stat, ks_p = ks_2samp(nonzero_fba, nonzero_pfba)
        print(f"KS D={ks_stat:.4f}, p={ks_p:.4f}")

    # -----------------------------------------------------------------------
    # 6. TOP DIFFERENCES
    # -----------------------------------------------------------------------
    print("\n6. REACTIONS WITH LARGEST DIFFERENCES")
    top_diff = df_fluxes.nlargest(10, "abs_difference")[["reaction", "FBA", "pFBA", "difference"]]
    print(top_diff.to_string(index=False))

    # -----------------------------------------------------------------------
    # PLOTTING
    # -----------------------------------------------------------------------
    fig, axes = plt.subplots(2, 3, figsize=(18, 11))
    fig.suptitle("FBA vs pFBA Comparison", fontsize=15)

    # Scatter
    axes[0, 0].scatter(df_fluxes["FBA"], df_fluxes["pFBA"], alpha=0.6, s=15)
    axes[0, 0].plot([-100, 100], [-100, 100], "r--", lw=1)
    axes[0, 0].set_xlabel("FBA flux")
    axes[0, 0].set_ylabel("pFBA flux")
    axes[0, 0].set_title(f"Flux correlation (r={pearson_r:.3f})")
    axes[0, 0].grid(alpha=0.3)

    # Histogram of differences
    axes[0, 1].hist(df_fluxes["difference"], bins=50, edgecolor="black")
    axes[0, 1].set_xlabel("FBA − pFBA")
    axes[0, 1].set_title("Flux difference distribution")
    axes[0, 1].axvline(0, color="red", ls="--")

    # Boxplot |flux|
    active_fba  = df_fluxes.loc[df_fluxes["active_FBA"],  "abs_FBA"]
    active_pfba = df_fluxes.loc[df_fluxes["active_pFBA"], "abs_pFBA"]
    axes[0, 2].boxplot([active_fba, active_pfba], labels=["FBA", "pFBA"])
    axes[0, 2].set_yscale("log")
    axes[0, 2].set_title("|Flux| (active reactions)")

    # Active counts
    axes[1, 0].bar(["FBA", "pFBA"], [active_fba.size, active_pfba.size],
                   color=["skyblue", "lightcoral"])
    for i, v in enumerate([active_fba.size, active_pfba.size]):
        axes[1, 0].text(i, v + 5, str(v), ha="center")
    axes[1, 0].set_ylabel("Active reactions")

    # CDF of |flux|
    sorted_fba  = np.sort(active_fba)
    sorted_pfba = np.sort(active_pfba)
    axes[1, 1].plot(sorted_fba,  np.linspace(0, 1, sorted_fba.size),  label="FBA")
    axes[1, 1].plot(sorted_pfba, np.linspace(0, 1, sorted_pfba.size), label="pFBA")
    axes[1, 1].set_xscale("log")
    axes[1, 1].set_xlabel("|Flux|")
    axes[1, 1].set_ylabel("Cumulative prob.")
    axes[1, 1].legend()

    # Heatmap of top 20 differences
    top_diff_reactions = df_fluxes.nlargest(20, 'abs_difference')
    heatmap_data = top_diff_reactions[['FBA', 'pFBA']].T
    vmax = np.abs(heatmap_data.values).max()
    im = axes[1,2].imshow(heatmap_data, aspect='auto', cmap='RdBu_r', 
                        vmin=-vmax, vmax=vmax)

    axes[1,2].set_title('Top 20 Most Different Reactions')
    axes[1,2].set_ylabel('Method')
    axes[1,2].set_xticks(np.arange(heatmap_data.shape[1]))
    axes[1,2].set_xticklabels(top_diff_reactions['reaction'], rotation=90, fontsize=8)
    axes[1,2].set_xlabel('Reaction')
    axes[1,2].set_yticks([0, 1])
    axes[1,2].set_yticklabels(['FBA', 'pFBA'])
    plt.colorbar(im, ax=axes[1,2], label='Flux Value')

    plt.tight_layout()
    plt.show()

    return df_fluxes

