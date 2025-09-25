# t-SNE Analysis: In Vivo vs In Vitro Enzyme-Substrate Space

## Overview

This pipeline creates a complete t-SNE visualization to compare in vivo and in vitro enzyme datasets, allowing you to analyze:

- **Coverage**: Are in vitro and in vivo entries populating the same biochemical neighborhoods?
- **Gaps**: Which regions have only in vivo (no assays) or only in vitro (no in vivo estimate)?
- **Local comparison**: In regions where both exist, do in vivo values run systematically lower than in vitro?

## Pipeline Components

### 1. Data Requirements

The pipeline expects three CSV files:
- `all_seqs_smiles.csv`: All unique (sequence, SMILES) pairs
- `kmax_clean.csv`: In vivo enzyme data with columns: sequence, SMILES, kcat_invivo
- `catpred_clean.csv`: In vitro enzyme data with columns: sequence, SMILES, kcat_invitro

### 2. Embedding Generation

**Protein Embeddings (1280D)**
- Uses ESM-2 650M model for protein sequence embeddings
- Automatically caches embeddings per unique sequence to avoid recomputation
- Mean pooling of residue-level embeddings

**Substrate Embeddings (2048D)**
- Uses RDKit ECFP6 (Morgan count fingerprint, radius=3)
- Handles invalid SMILES gracefully with zero vectors
- Count-based fingerprints for better chemical similarity representation

### 3. Dimensionality Reduction

**Feature Combination & Scaling**
- Concatenates protein (1280D) + substrate (2048D) = 3328D
- StandardScaler normalization for stable training

**PCA Reduction**
- Reduces to 50-128 dimensions (configurable, default: 100)
- Denoises data and speeds up t-SNE computation
- Reports explained variance ratio

**t-SNE Mapping**
- Maps to 2D for visualization
- Uses optimized parameters (perplexity=30, n_iter=1000)
- Consistent random seed for reproducibility

### 4. Visualization & Analysis

**Main t-SNE Plot**
- Red points: In vitro data (CatPred)
- Blue points: In vivo data (kmax)
- Clusters indicate similar enzyme-substrate combinations

**Coverage Analysis**
- Grid-based spatial analysis of the t-SNE space
- Identifies regions with:
  - Only in vitro data
  - Only in vivo data
  - Both datasets
  - No coverage
- Quantifies overlap and gap statistics

## Installation

### Required Dependencies

```bash
# Core dependencies (already in requirements.txt)
pip install numpy pandas scikit-learn matplotlib seaborn

# Machine learning models
pip install fair-esm torch

# Chemical informatics
pip install rdkit-pypi

# Optional: Interactive visualization
pip install plotly
```

### Quick Install
```bash
pip install -r requirements.txt
```

## Usage

### Option 1: Jupyter Notebook (Recommended)

1. Open `experiments/tsne_invivo_invitro.ipynb`
2. Run all cells sequentially
3. Results automatically saved to `results/tsne_analysis/`

### Option 2: Direct Script Execution

```python
from scripts.tsne_kcats import run_tsne_pipeline

results = run_tsne_pipeline(
    all_seqs_file="data/processed/experiments/all_seqs_smiles.csv",
    kmax_file="data/processed/experiments/kmax_clean.csv",
    catpred_file="data/processed/experiments/catpred_clean.csv",
    output_dir="results/tsne_analysis",
    pca_components=100,
    cache_dir="cache"
)
```

### Option 3: Command Line

```bash
cd scripts
python tsne_kcats.py
```

## Testing

Run the test suite to verify everything is working:

```bash
python test_tsne_pipeline.py
```

This will check:
- All required packages are installed
- Data files exist
- Pipeline functions import correctly
- Basic embedding computation works

## Output Files

The pipeline generates:

```
results/tsne_analysis/
├── tsne_invivo_vs_invitro.png     # Main t-SNE visualization
├── coverage_analysis.png          # Coverage analysis plots
├── tsne_results.pkl               # Complete results object
└── analysis_summary.txt           # Text summary of findings
```

## Interpreting Results

### t-SNE Plot Analysis

**Clustering Patterns**
- Tight clusters = similar enzyme-substrate combinations
- Isolated points = unique biochemical space
- Color mixing = regions where both datasets have data

**Dataset Distribution**
- Red regions = only in vitro measurements available
- Blue regions = only in vivo estimates available
- Purple regions = both datasets present

### Coverage Statistics

**High Overlap (>30%)**
- Datasets cover similar biochemical space
- Good agreement on important enzyme-substrate pairs
- Focus analysis on systematic differences in overlapping regions

**Low Overlap (<10%)**
- Datasets occupy different biochemical niches
- Large gaps in coverage for both datasets
- Priority: expand coverage in complementary directions

**Gap Analysis**
- In vitro only regions → need in vivo measurements
- In vivo only regions → need in vitro assays
- No coverage regions → unexplored biochemical space

### Research Questions Answered

1. **Coverage**: Overlap percentage indicates shared biochemical neighborhoods
2. **Gaps**: Visual and quantitative identification of single-dataset regions
3. **Local comparison**: Framework established for comparing kcat values in overlapping regions

## Advanced Usage

### Parameter Tuning

**PCA Components**
```python
# More components = retain more information, slower t-SNE
# Fewer components = faster, may lose important features
results = run_tsne_pipeline(..., pca_components=128)  # vs default 100
```

**t-SNE Parameters**
```python
# Modify in create_tsne_plot function:
# - perplexity: 5-50 (smaller = more local structure)
# - n_iter: 1000+ (more iterations = better convergence)
```

**Grid Resolution**
```python
# Modify in analyze_coverage function:
# - grid_size: 50+ (higher = finer spatial resolution)
```

### Custom Analysis

**Enzyme Class Mapping**
```python
# Add EC number or enzyme family information
# Color code by functional class instead of dataset
```

**Substrate Chemical Properties**
```python
# Include additional molecular descriptors
# Analyze chemical diversity across regions
```

**Statistical Comparison**
```python
# In overlapping regions, compare kcat distributions
# Test for systematic in vivo vs in vitro differences
```

## Performance Considerations

### Computational Requirements

**Memory Usage**
- ESM-2 model: ~3GB GPU/CPU memory
- Embeddings: ~50MB per 1000 sequences
- Full pipeline: 4-8GB RAM recommended

**Processing Time**
- Protein embeddings: ~1-2 sec per sequence (cached after first run)
- Substrate embeddings: ~0.1 sec per SMILES
- PCA + t-SNE: 1-5 minutes for 1000-2000 samples

### Optimization Tips

1. **Use Caching**: Protein embeddings are automatically cached
2. **Batch Processing**: Pipeline handles full datasets efficiently
3. **GPU Acceleration**: ESM-2 automatically uses GPU if available
4. **PCA First**: Always apply PCA before t-SNE for speed

## Troubleshooting

### Common Issues

**Import Errors**
```bash
# ESM not found
pip install fair-esm

# RDKit not found  
pip install rdkit-pypi

# CUDA issues
pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cpu
```

**Memory Issues**
```python
# Reduce PCA components
results = run_tsne_pipeline(..., pca_components=50)

# Process in smaller batches
# Modify get_protein_embeddings to process sequences in batches
```

**Invalid SMILES**
- Pipeline handles invalid SMILES gracefully
- Check data quality if many warnings appear
- Consider pre-filtering or cleaning SMILES strings

### Performance Issues

**Slow Protein Embedding**
- First run is slow (downloads ESM-2 model)
- Subsequent runs use cached embeddings
- Consider using smaller ESM model for testing

**t-SNE Convergence**
- Increase n_iter if plot looks unstructured
- Try different perplexity values
- Check PCA explained variance (should be >80%)

## Extensions and Future Work

### Immediate Extensions

1. **UMAP Alternative**: Add UMAP dimensionality reduction for comparison
2. **Interactive Plots**: Use plotly for explorable visualizations
3. **Quantitative Comparison**: Statistical tests for in vivo vs in vitro differences
4. **Cluster Analysis**: Identify and characterize distinct biochemical regions

### Research Applications

1. **Enzyme Engineering**: Identify underexplored enzyme-substrate combinations
2. **Assay Design**: Prioritize in vitro assays for in vivo-only regions
3. **Model Validation**: Compare prediction accuracy across different regions
4. **Database Curation**: Identify systematic biases in data collection

### Integration with Other Tools

1. **GEM Integration**: Map results to metabolic network context
2. **Phylogenetic Analysis**: Correlate with evolutionary relationships
3. **Structural Analysis**: Include protein structure features
4. **Kinetic Modeling**: Use insights for parameter estimation

---

## Support

For questions or issues:
1. Check the test suite: `python test_tsne_pipeline.py`
2. Review the Jupyter notebook examples
3. Examine output logs for specific error messages
4. Consider the troubleshooting section above

The pipeline is designed to be robust and provide clear error messages to guide troubleshooting.

