"""
Module Description: tsne_kcats.py

Purpose: 
Complete pipeline for comparing in vivo vs in vitro enzyme datasets using t-SNE visualization.
Creates embeddings from protein sequences (ESM-2) and substrate SMILES (ECFP6), then visualizes
the enzyme-substrate space to identify coverage gaps and local differences.

Functions:
- get_protein_embeddings: Get ESM-2 650M protein embeddings with caching
- get_substrate_embeddings: Get RDKit ECFP6 substrate embeddings  
- create_combined_embeddings: Concatenate and scale protein + substrate embeddings
- apply_pca_reduction: Apply PCA dimensionality reduction for denoising
- create_tsne_plot: Create t-SNE visualization colored by dataset
- analyze_coverage: Analyze coverage, gaps, and local comparisons
- run_tsne_pipeline: Complete pipeline execution
"""

import os
import pickle
import hashlib
import numpy as np
import pandas as pd
from typing import Dict, List, Tuple, Optional
import warnings
import matplotlib.lines as mlines
warnings.filterwarnings('ignore')

# ML and visualization libraries
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.manifold import TSNE
import matplotlib.pyplot as plt
import seaborn as sns

# Chemical informatics
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
import torch

# ESM for protein embeddings
try:
    import esm
    ESM_AVAILABLE = True
except ImportError:
    ESM_AVAILABLE = False
    print("ESM not available. Install with: pip install fair-esm")


def get_protein_embeddings(sequences: List[str], cache_dir: str = "cache/protein_embeddings") -> np.ndarray:
    """
    Get ESM-2 650M protein embeddings with caching to avoid recomputation.
    """
    if not ESM_AVAILABLE:
        raise ImportError("ESM library not available. Install with: pip install fair-esm")
    
    # Create cache directory
    os.makedirs(cache_dir, exist_ok=True)
    
    # Load ESM-2 650M model
    model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
    batch_converter = alphabet.get_batch_converter()
    model.eval()
    
    embeddings = []
    
    for seq in sequences:
        # Create hash for sequence to use as cache key
        seq_hash = hashlib.md5(seq.encode()).hexdigest()
        cache_file = os.path.join(cache_dir, f"{seq_hash}.pkl")
        
        # Try to load from cache
        if os.path.exists(cache_file):
            print(f"Loading embedding from cache: {cache_file}")
            with open(cache_file, 'rb') as f:
                embedding = pickle.load(f)
        else:
            print(f"No cache found. Computing embedding for sequence: {seq}")
            # Compute embedding
            batch_labels, batch_strs, batch_tokens = batch_converter([("protein", seq)])
            
            with torch.no_grad():
                results = model(batch_tokens, repr_layers=[33], return_contacts=False)
            
            # Extract per-sequence embedding (mean pooling of residue embeddings)
            token_embeddings = results["representations"][33]
            # Remove batch dimension and exclude special tokens
            sequence_embedding = token_embeddings[0, 1 : len(seq) + 1].mean(0)
            embedding = sequence_embedding.numpy()
            
            # Cache the embedding
            with open(cache_file, 'wb') as f:
                pickle.dump(embedding, f)
        
        embeddings.append(embedding)
    
    return np.array(embeddings)


def get_substrate_embeddings(smiles: List[str], radius: int = 3, n_bits: int = 2048) -> np.ndarray:
    """
    Get RDKit ECFP6 (Morgan fingerprint) substrate embeddings.
    """
    embeddings = []
    
    # Suppress RDKit warnings for cleaner output
    from rdkit import RDLogger
    RDLogger.DisableLog('rdApp.*')
    
    # MorganGenerator API if available, fallback to old API
    try:
        from rdkit.Chem.rdMolDescriptors import MorganGenerator
        use_new_api = True
        mg = MorganGenerator(radius=radius, fpSize=n_bits)
        print("Using new MorganGenerator API (no deprecation warnings)")
    except (ImportError, AttributeError):
        use_new_api = False
        print("Using legacy Morgan fingerprint API")
    
    for smi in smiles:
        try:
            mol = Chem.MolFromSmiles(smi)
            if mol is not None:
                if use_new_api:
                    # Use new MorganGenerator API
                    fingerprint = mg.GetFingerprint(mol)
                    # Convert to numpy array
                    arr = np.zeros((n_bits,), dtype=np.float32)
                    for i in range(n_bits):
                        arr[i] = fingerprint[i]
                    embeddings.append(arr)
                else:
                    # Fallback to old API
                    fingerprint = rdMolDescriptors.GetMorganFingerprintAsBitVect(
                        mol, radius=radius, nBits=n_bits
                    )
                    # Convert to numpy array
                    arr = np.zeros((n_bits,), dtype=np.float32)
                    for i in range(n_bits):
                        arr[i] = fingerprint[i]
                    embeddings.append(arr)
            else:
                # Invalid SMILES - use zero vector
                embeddings.append(np.zeros(n_bits, dtype=np.float32))
                print(f"Warning: Invalid SMILES: {smi}")
        except Exception as e:
            print(f"Error processing SMILES {smi}: {e}")
            embeddings.append(np.zeros(n_bits, dtype=np.float32))
    
    return np.array(embeddings)


def create_combined_embeddings(protein_embeddings: np.ndarray, 
                             substrate_embeddings: np.ndarray) -> Tuple[np.ndarray, StandardScaler]:
    """
    Concatenate and scale protein + substrate embeddings.
    """
    # Concatenate embeddings
    combined = np.hstack([protein_embeddings, substrate_embeddings])
    
    # Scale features
    scaler = StandardScaler()
    combined_scaled = scaler.fit_transform(combined)
    
    return combined_scaled, scaler


def apply_pca_reduction(embeddings: np.ndarray, n_components: int = 100) -> Tuple[np.ndarray, PCA]:
    """
    Apply PCA dimensionality reduction for denoising and speedup.
    """
    pca = PCA(n_components=n_components, random_state=42)
    embeddings_pca = pca.fit_transform(embeddings)
    
    print(f"PCA explained variance ratio: {pca.explained_variance_ratio_.sum():.3f}")
    
    return embeddings_pca, pca


def create_tsne_plot(embeddings: np.ndarray, 
                    labels: np.ndarray,
                    dataset_names: List[str],
                    title: str = "t-SNE: In Vivo vs In Vitro Enzyme-Substrate Space",
                    figsize: Tuple[int, int] = (12, 8),
                    save_path: Optional[str] = None,
                    subsystem_info: Optional[np.ndarray] = None) -> plt.Figure:
    """
    Create t-SNE visualization colored by dataset and subsystem.
    """
    # Run t-SNE
    print("Running t-SNE...")
    # Handle different scikit-learn versions (n_iter vs max_iter)
    try:
        tsne = TSNE(n_components=2, random_state=42, perplexity=30, max_iter=1000)
    except TypeError:
        # Fallback for older scikit-learn versions
        tsne = TSNE(n_components=2, random_state=42, perplexity=30, n_iter=1000)
    
    embeddings_2d = tsne.fit_transform(embeddings)
    
    # Create plot
    fig, ax = plt.subplots(figsize=figsize)
    
    # Plot in vitro data (crosses)
    in_vitro_mask = labels == 0
    if np.any(in_vitro_mask):
        ax.scatter(embeddings_2d[in_vitro_mask, 0], embeddings_2d[in_vitro_mask, 1], 
                  c='#808080', marker='x', label=dataset_names[0], alpha=0.6, s=30, linewidths=2)
    
    # Plot in vivo data (circles, colored by subsystem)
    in_vivo_mask = labels == 1
    if np.any(in_vivo_mask) and subsystem_info is not None:
        # Get unique subsystems for in vivo data
        vivo_subsystems = subsystem_info[in_vivo_mask]
        unique_subsystems = np.unique(vivo_subsystems)
        
        # Create color palette for subsystems
        from matplotlib.colors import ListedColormap
        import matplotlib.cm as cm
        cmap = cm.get_cmap('tab20')  # Good for categorical data
        subsystem_colors = [cmap(i / len(unique_subsystems)) for i in range(len(unique_subsystems))]
        
        # Plot each subsystem
        for i, subsystem in enumerate(unique_subsystems):
            subsystem_mask = (labels == 1) & (subsystem_info == subsystem)
            if np.any(subsystem_mask):
                ax.scatter(embeddings_2d[subsystem_mask, 0], embeddings_2d[subsystem_mask, 1], 
                          c=[subsystem_colors[i]], marker='o', label=f'{subsystem}', 
                          alpha=0.7, s=25, linewidths=0.5)
    
    elif np.any(in_vivo_mask):
        # Fallback: plot all in vivo as blue circles if no subsystem info
        ax.scatter(embeddings_2d[in_vivo_mask, 0], embeddings_2d[in_vivo_mask, 1], 
                  c='#3498DB', marker='o', label=dataset_names[1], alpha=0.6, s=25)
    
    ax.set_xlabel('t-SNE 1')
    ax.set_ylabel('t-SNE 2')
    ax.set_title(title, fontstyle='italic')
    
    # Create legend with manageable number of entries
    handles, labels_legend = ax.get_legend_handles_labels()

    # Create an empty (invisible) handle for the extra line
    dummy_handle = mlines.Line2D([], [], color='none')

    # Insert the new label right after 'In vitro'
    if labels_legend and labels_legend[0] == 'In vitro':
        handles.insert(1, dummy_handle)
        labels_legend.insert(1, 'In vivo:')

    # Build the legend
    if len(handles) > 15:
        ax.legend(handles[:15], labels_legend[:15],
                bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.figtext(0.99, 0.01,
                    f"Showing first 15 of {len(handles)} categories",
                    ha='right', va='bottom', fontsize=8)
    else:
        ax.legend(handles, labels_legend, loc='lower right')

    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Plot saved to: {save_path}")
    
    return fig


def analyze_coverage(embeddings_2d: np.ndarray, 
                    labels: np.ndarray,
                    dataset_names: List[str],
                    grid_size: int = 50) -> Dict:
    """
    Analyze coverage, gaps, and local comparisons in the t-SNE space.
    """
    # Create grid
    x_min, x_max = embeddings_2d[:, 0].min(), embeddings_2d[:, 0].max()
    y_min, y_max = embeddings_2d[:, 1].min(), embeddings_2d[:, 1].max()
    
    x_bins = np.linspace(x_min, x_max, grid_size)
    y_bins = np.linspace(y_min, y_max, grid_size)
    
    # Count points in each grid cell by dataset
    coverage_maps = {}
    for i, dataset_name in enumerate(dataset_names):
        mask = labels == i
        points = embeddings_2d[mask]
        coverage_map, _, _ = np.histogram2d(points[:, 0], points[:, 1], 
                                          bins=[x_bins, y_bins])
        coverage_maps[dataset_name] = coverage_map
    
    # Calculate coverage statistics
    in_vitro_coverage = coverage_maps[dataset_names[0]] > 0
    in_vivo_coverage = coverage_maps[dataset_names[1]] > 0
    
    # Areas
    total_area = grid_size * grid_size
    in_vitro_only = in_vitro_coverage & ~in_vivo_coverage
    in_vivo_only = in_vivo_coverage & ~in_vitro_coverage
    both_coverage = in_vitro_coverage & in_vivo_coverage
    no_coverage = ~(in_vitro_coverage | in_vivo_coverage)
    
    results = {
        'total_cells': total_area,
        'in_vitro_only_cells': in_vitro_only.sum(),
        'in_vivo_only_cells': in_vivo_only.sum(),
        'both_datasets_cells': both_coverage.sum(),
        'no_coverage_cells': no_coverage.sum(),
        'in_vitro_only_fraction': in_vitro_only.sum() / total_area,
        'in_vivo_only_fraction': in_vivo_only.sum() / total_area,
        'both_datasets_fraction': both_coverage.sum() / total_area,
        'no_coverage_fraction': no_coverage.sum() / total_area,
        'coverage_maps': coverage_maps,
        'grid_bounds': (x_bins, y_bins)
    }
    
    return results


def analyze_subsystem_coverage(embeddings_2d: np.ndarray, 
                              labels: np.ndarray,
                              subsystems: np.ndarray,
                              dataset_names: List[str],
                              grid_size: int = 50) -> Dict:
    """
    Analyze subsystem proportions in different coverage regions.
    """
    # Get basic coverage analysis
    coverage_results = analyze_coverage(embeddings_2d, labels, dataset_names, grid_size)
    
    # Use the same grid bounds as the coverage analysis
    x_bins, y_bins = coverage_results['grid_bounds']
    
    # Get coverage maps
    in_vitro_coverage = coverage_results['coverage_maps'][dataset_names[0]] > 0
    in_vivo_coverage = coverage_results['coverage_maps'][dataset_names[1]] > 0
    
    # Define regions
    in_vivo_only_region = in_vivo_coverage & ~in_vitro_coverage
    both_coverage_region = in_vitro_coverage & in_vivo_coverage
    
    # Get in vivo points and their subsystems
    in_vivo_mask = labels == 1
    in_vivo_points = embeddings_2d[in_vivo_mask]
    in_vivo_subsystems = subsystems[in_vivo_mask]
    
    # Map in vivo points to grid cells using numpy's histogram2d method
    # This ensures exact consistency with how coverage maps were created
    in_vivo_counts, _, _ = np.histogram2d(in_vivo_points[:, 0], in_vivo_points[:, 1], 
                                         bins=[x_bins, y_bins])
    
    # For each in vivo point, find which grid cell it belongs to
    subsystems_in_vivo_only = []
    subsystems_in_both = []
    
    for i, point in enumerate(in_vivo_points):
        # Find grid indices for this point
        x_idx = np.searchsorted(x_bins[1:], point[0])
        y_idx = np.searchsorted(y_bins[1:], point[1])
        
        # Ensure indices are within bounds
        x_idx = min(x_idx, grid_size - 1)
        y_idx = min(y_idx, grid_size - 1)
        
        # Classify by region
        if in_vivo_only_region[x_idx, y_idx]:
            subsystems_in_vivo_only.append(in_vivo_subsystems[i])
        elif both_coverage_region[x_idx, y_idx]:
            subsystems_in_both.append(in_vivo_subsystems[i])
    
    # Count subsystem proportions
    from collections import Counter
    
    vivo_only_counts = Counter(subsystems_in_vivo_only)
    both_counts = Counter(subsystems_in_both)
    
    return {
        'coverage_results': coverage_results,
        'vivo_only_subsystems': dict(vivo_only_counts),
        'both_coverage_subsystems': dict(both_counts),
        'vivo_only_total': len(subsystems_in_vivo_only),
        'both_coverage_total': len(subsystems_in_both)
    }


def create_subsystem_pie_charts(subsystem_analysis: Dict,
                               figsize: Tuple[int, int] = (15, 6),
                               save_path: Optional[str] = None) -> plt.Figure:
    """
    Create horizontal plot with two pie charts showing subsystem proportions.
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)
    
    # Color palette for consistency
    import matplotlib.cm as cm
    all_subsystems = set(subsystem_analysis['vivo_only_subsystems'].keys()) | \
                    set(subsystem_analysis['both_coverage_subsystems'].keys())
    cmap = cm.get_cmap('tab20')
    subsystem_colors = {subsystem: cmap(i / len(all_subsystems)) 
                       for i, subsystem in enumerate(sorted(all_subsystems))}
    
    # Pie chart 1: In vivo only regions
    vivo_only_data = subsystem_analysis['vivo_only_subsystems']
    if vivo_only_data:
        labels1 = list(vivo_only_data.keys())
        sizes1 = list(vivo_only_data.values())
        colors1 = [subsystem_colors[label] for label in labels1]
        
        wedges1, texts1, autotexts1 = ax1.pie(sizes1, labels=labels1, colors=colors1, 
                                              autopct='%1.1f%%', startangle=90)
        
        # Improve text readability
        for autotext in autotexts1:
            autotext.set_color('white')
            autotext.set_fontweight('bold')
            autotext.set_fontsize(8)
        
        # Rotate labels for better readability
        for text in texts1:
            text.set_fontsize(8)
            text.set_rotation_mode('anchor')
    else:
        ax1.text(0.5, 0.5, 'No data', ha='center', va='center', transform=ax1.transAxes)
    
    ax1.set_title(f'In Vivo Only Regions\n({subsystem_analysis["vivo_only_total"]} enzymes)', 
                  fontsize=12, fontweight='bold')
    
    # Pie chart 2: Both coverage regions
    both_data = subsystem_analysis['both_coverage_subsystems']
    if both_data:
        labels2 = list(both_data.keys())
        sizes2 = list(both_data.values())
        colors2 = [subsystem_colors[label] for label in labels2]
        
        wedges2, texts2, autotexts2 = ax2.pie(sizes2, labels=labels2, colors=colors2,
                                              autopct='%1.1f%%', startangle=90)
        
        # Improve text readability
        for autotext in autotexts2:
            autotext.set_color('white')
            autotext.set_fontweight('bold')
            autotext.set_fontsize(8)
        
        # Rotate labels for better readability
        for text in texts2:
            text.set_fontsize(8)
            text.set_rotation_mode('anchor')
    else:
        ax2.text(0.5, 0.5, 'No data', ha='center', va='center', transform=ax2.transAxes)
    
    ax2.set_title(f'Both Datasets Regions\n({subsystem_analysis["both_coverage_total"]} enzymes)', 
                  fontsize=12, fontweight='bold')
    
    # Overall title
    fig.suptitle('Subsystem Distribution by Coverage Region', fontsize=14, fontweight='bold')
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Subsystem pie charts saved to: {save_path}")
    
    return fig


def create_coverage_plot(coverage_results: Dict, 
                        dataset_names: List[str],
                        save_path: Optional[str] = None) -> plt.Figure:
    """
    Create coverage analysis visualization.
    """
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    
    x_bins, y_bins = coverage_results['grid_bounds']
    
    # Plot 1: In vitro coverage
    im1 = axes[0, 0].imshow(coverage_results['coverage_maps'][dataset_names[0]].T, 
                           origin='lower', extent=[x_bins[0], x_bins[-1], y_bins[0], y_bins[-1]],
                           cmap='Reds', alpha=0.8)
    axes[0, 0].set_title(f'{dataset_names[0]} Coverage')
    axes[0, 0].set_xlabel('t-SNE 1')
    axes[0, 0].set_ylabel('t-SNE 2')
    plt.colorbar(im1, ax=axes[0, 0])
    
    # Plot 2: In vivo coverage
    im2 = axes[0, 1].imshow(coverage_results['coverage_maps'][dataset_names[1]].T, 
                           origin='lower', extent=[x_bins[0], x_bins[-1], y_bins[0], y_bins[-1]],
                           cmap='Blues', alpha=0.8)
    axes[0, 1].set_title(f'{dataset_names[1]} Coverage')
    axes[0, 1].set_xlabel('t-SNE 1')
    axes[0, 1].set_ylabel('t-SNE 2')
    plt.colorbar(im2, ax=axes[0, 1])
    
    # Plot 3: Coverage overlap
    in_vitro_coverage = coverage_results['coverage_maps'][dataset_names[0]] > 0
    in_vivo_coverage = coverage_results['coverage_maps'][dataset_names[1]] > 0
    
    overlap_map = np.zeros_like(in_vitro_coverage, dtype=int)
    overlap_map[in_vitro_coverage & ~in_vivo_coverage] = 1  # In vitro only
    overlap_map[in_vivo_coverage & ~in_vitro_coverage] = 2   # In vivo only
    overlap_map[in_vitro_coverage & in_vivo_coverage] = 3    # Both
    
    colors = ['white', '#E74C3C', '#3498DB', '#9B59B6']  # White, Red, Blue, Purple
    from matplotlib.colors import ListedColormap
    cmap = ListedColormap(colors)
    
    im3 = axes[1, 0].imshow(overlap_map.T, origin='lower', 
                           extent=[x_bins[0], x_bins[-1], y_bins[0], y_bins[-1]],
                           cmap=cmap, vmin=0, vmax=3)
    axes[1, 0].set_title('Coverage Overlap')
    axes[1, 0].set_xlabel('t-SNE 1')
    axes[1, 0].set_ylabel('t-SNE 2')
    
    # Custom colorbar for overlap plot
    cbar = plt.colorbar(im3, ax=axes[1, 0], ticks=[0.375, 1.125, 1.875, 2.625])
    cbar.ax.set_yticklabels(['None', 'In vitro only', 'In vivo only', 'Both'])
    
    # Plot 4: Coverage statistics
    categories = ['In vitro only', 'In vivo only', 'Both datasets', 'No coverage']
    fractions = [coverage_results['in_vitro_only_fraction'],
                coverage_results['in_vivo_only_fraction'], 
                coverage_results['both_datasets_fraction'],
                coverage_results['no_coverage_fraction']]
    
    bars = axes[1, 1].bar(categories, fractions, 
                         color=['#E74C3C', '#3498DB', '#9B59B6', '#95A5A6'])
    axes[1, 1].set_title('Coverage Statistics')
    axes[1, 1].set_ylabel('Fraction of Space')
    axes[1, 1].tick_params(axis='x', rotation=45)
    
    # Add value labels on bars
    for bar, fraction in zip(bars, fractions):
        height = bar.get_height()
        axes[1, 1].text(bar.get_x() + bar.get_width()/2., height + 0.01,
                        f'{fraction:.3f}', ha='center', va='bottom')
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Coverage plot saved to: {save_path}")
    
    return fig


def run_tsne_pipeline(all_seqs_file: str,
                     kmax_file: str, 
                     catpred_file: str,
                     output_dir: str = "results/tsne_analysis",
                     pca_components: int = 100,
                     cache_dir: str = "cache",
                     use_subsystems: bool = True) -> Dict:
    """
    Complete t-SNE pipeline execution.
    """
    print("=== Starting t-SNE Pipeline ===")
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Load data
    print("Loading datasets...")
    all_seqs = pd.read_csv(all_seqs_file)
    kmax_data = pd.read_csv(kmax_file)  # in vivo
    catpred_data = pd.read_csv(catpred_file)  # in vitro
    
    print(f"All sequences: {len(all_seqs)} entries")
    print(f"In vivo (kmax): {len(kmax_data)} entries")
    print(f"In vitro (catpred): {len(catpred_data)} entries")
    
    # Step 1: Get protein embeddings
    print("\n=== Step 1: Computing protein embeddings ===")
    unique_sequences = all_seqs['sequence'].unique()
    print(f"Computing embeddings for {len(unique_sequences)} unique sequences...")
    
    protein_embeddings = get_protein_embeddings(unique_sequences, 
                                               cache_dir=os.path.join(cache_dir, "protein_embeddings"))
    
    # Create mapping from sequence to embedding
    seq_to_embedding = dict(zip(unique_sequences, protein_embeddings))
    
    # Map embeddings to all sequences
    all_protein_embeddings = np.array([seq_to_embedding[seq] for seq in all_seqs['sequence']])
    
    # Step 2: Get substrate embeddings
    print("\n=== Step 2: Computing substrate embeddings ===")
    substrate_embeddings = get_substrate_embeddings(all_seqs['SMILES'].tolist())
    
    # Step 3: Combine and scale embeddings
    print("\n=== Step 3: Combining and scaling embeddings ===")
    combined_embeddings, scaler = create_combined_embeddings(all_protein_embeddings, substrate_embeddings)
    
    print(f"Combined embeddings shape: {combined_embeddings.shape}")
    
    # Step 4: PCA reduction
    print("\n=== Step 4: Applying PCA reduction ===")
    pca_embeddings, pca = apply_pca_reduction(combined_embeddings, n_components=pca_components)
    
    # Step 5: Create labels for datasets and subsystem information
    print("\n=== Step 5: Creating dataset labels and subsystem mapping ===")
    
    # Create mappings to identify which dataset each entry belongs to
    kmax_key_to_subsystem = {}
    kmax_key_set = set()
    for _, row in kmax_data.iterrows():
        key = (row['sequence'], row['SMILES'])
        kmax_key_set.add(key)
        if use_subsystems and 'subsystem' in row:
            kmax_key_to_subsystem[key] = row['subsystem']
    
    catpred_key_set = set()
    for _, row in catpred_data.iterrows():
        catpred_key_set.add((row['sequence'], row['SMILES']))
    
    # Create labels: 0 = in vitro, 1 = in vivo
    labels = []
    subsystems = []
    for _, row in all_seqs.iterrows():
        key = (row['sequence'], row['SMILES'])
        in_kmax = key in kmax_key_set
        in_catpred = key in catpred_key_set
        
        if in_catpred and not in_kmax:
            labels.append(0)  # In vitro only
            subsystems.append(None)  # No subsystem for in vitro
        elif in_kmax and not in_catpred:
            labels.append(1)  # In vivo only
            subsystems.append(kmax_key_to_subsystem.get(key, 'Unknown'))
        elif in_kmax and in_catpred:
            labels.append(1)  # Both - assign to in vivo for visualization
            subsystems.append(kmax_key_to_subsystem.get(key, 'Unknown'))
        else:
            labels.append(-1)  # Neither (shouldn't happen with proper data)
            subsystems.append(None)
    
    labels = np.array(labels)
    subsystems = np.array(subsystems)
    
    # Filter out entries that don't belong to either dataset
    valid_mask = labels >= 0
    pca_embeddings_filtered = pca_embeddings[valid_mask]
    labels_filtered = labels[valid_mask]
    subsystems_filtered = subsystems[valid_mask]
    
    print(f"In vitro entries: {np.sum(labels_filtered == 0)}")
    print(f"In vivo entries: {np.sum(labels_filtered == 1)}")
    
    if use_subsystems:
        unique_subsystems = np.unique(subsystems_filtered[subsystems_filtered != None])
        print(f"Unique subsystems in in vivo data: {len(unique_subsystems)}")
        print(f"Subsystems: {list(unique_subsystems)}")
    
    # Step 6: Create t-SNE plot
    print("\n=== Step 6: Creating t-SNE visualization ===")
    dataset_names = ['In vitro', 'In vivo']
    
    tsne_fig = create_tsne_plot(pca_embeddings_filtered, 
                               labels_filtered,
                               dataset_names,
                               save_path=os.path.join(output_dir, "tsne_invivo_vs_invitro_subsystems.png"),
                               subsystem_info=subsystems_filtered if use_subsystems else None)
    
    # Get 2D embeddings for analysis
    # Handle different scikit-learn versions (n_iter vs max_iter)
    try:
        tsne = TSNE(n_components=2, random_state=42, perplexity=30, max_iter=1000)
    except TypeError:
        # Fallback for older scikit-learn versions
        tsne = TSNE(n_components=2, random_state=42, perplexity=30, n_iter=1000)
    
    embeddings_2d = tsne.fit_transform(pca_embeddings_filtered)
    
    # Step 7: Coverage analysis
    print("\n=== Step 7: Analyzing coverage ===")
    coverage_results = analyze_coverage(embeddings_2d, labels_filtered, dataset_names)
    
    # Print coverage statistics
    print(f"Coverage Statistics:")
    print(f"  In vitro only: {coverage_results['in_vitro_only_fraction']:.1%}")
    print(f"  In vivo only: {coverage_results['in_vivo_only_fraction']:.1%}")  
    print(f"  Both datasets: {coverage_results['both_datasets_fraction']:.1%}")
    print(f"  No coverage: {coverage_results['no_coverage_fraction']:.1%}")
    
    # Create coverage plot
    coverage_fig = create_coverage_plot(coverage_results, dataset_names,
                                       save_path=os.path.join(output_dir, "coverage_analysis.png"))
    
    # Step 8: Subsystem analysis (if subsystems are available)
    subsystem_pie_fig = None
    if use_subsystems and subsystems_filtered is not None:
        print("\n=== Step 8: Analyzing subsystem distribution by coverage regions ===")
        subsystem_analysis = analyze_subsystem_coverage(embeddings_2d, labels_filtered, 
                                                       subsystems_filtered, dataset_names)
        
        print(f"Subsystem Statistics:")
        print(f"  In vivo only regions: {subsystem_analysis['vivo_only_total']} enzymes")
        print(f"  Both coverage regions: {subsystem_analysis['both_coverage_total']} enzymes")
        
        # Show top subsystems in each region
        vivo_only_subs = subsystem_analysis['vivo_only_subsystems']
        both_subs = subsystem_analysis['both_coverage_subsystems']
        
        if vivo_only_subs:
            print(f"  Top subsystems in vivo-only regions:")
            sorted_vivo = sorted(vivo_only_subs.items(), key=lambda x: x[1], reverse=True)
            for subsystem, count in sorted_vivo[:5]:
                print(f"    - {subsystem}: {count} enzymes ({count/subsystem_analysis['vivo_only_total']*100:.1f}%)")
        
        if both_subs:
            print(f"  Top subsystems in both-coverage regions:")
            sorted_both = sorted(both_subs.items(), key=lambda x: x[1], reverse=True)
            for subsystem, count in sorted_both[:5]:
                print(f"    - {subsystem}: {count} enzymes ({count/subsystem_analysis['both_coverage_total']*100:.1f}%)")
        
        # Create subsystem pie charts
        subsystem_pie_fig = create_subsystem_pie_charts(
            subsystem_analysis,
            save_path=os.path.join(output_dir, "subsystem_pie_charts.png")
        )
    
    # Save results
    results = {
        'all_seqs': all_seqs,
        'kmax_data': kmax_data,
        'catpred_data': catpred_data,
        'combined_embeddings': combined_embeddings,
        'pca_embeddings': pca_embeddings_filtered,
        'embeddings_2d': embeddings_2d,
        'labels': labels_filtered,
        'subsystems': subsystems_filtered if use_subsystems else None,
        'dataset_names': dataset_names,
        'coverage_results': coverage_results,
        'subsystem_analysis': subsystem_analysis if use_subsystems and subsystems_filtered is not None else None,
        'scaler': scaler,
        'pca': pca,
        'tsne_fig': tsne_fig,
        'coverage_fig': coverage_fig,
        'subsystem_pie_fig': subsystem_pie_fig
    }
    
    # Save embeddings and results
    results_file = os.path.join(output_dir, "tsne_results.pkl")
    with open(results_file, 'wb') as f:
        # Don't save the figures to avoid issues
        save_results = {k: v for k, v in results.items() if 'fig' not in k}
        pickle.dump(save_results, f)
    
    print(f"\n=== Pipeline Complete ===")
    print(f"Results saved to: {output_dir}")
    print(f"Key findings:")
    print(f"  - Total enzyme-substrate pairs analyzed: {len(labels_filtered)}")
    print(f"  - Enzyme-substrate space regions with both datasets: {coverage_results['both_datasets_fraction']:.1%}")
    print(f"  - Regions with only in vitro data: {coverage_results['in_vitro_only_fraction']:.1%}")
    print(f"  - Regions with only in vivo data: {coverage_results['in_vivo_only_fraction']:.1%}")
    
    return results


