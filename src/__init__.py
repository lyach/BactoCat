"""
BactoCat: A pipeline for building kapp datasets from genome-scale metabolic models.

This package provides tools for:
- Enzyme classification from GPR rules
- Gene sequence mapping via UniProt
- Substrate mapping to SMILES structures
- Flux analysis and kapp calculations

Modules
-------
config
    Configuration management and path definitions
enzyme_classifier
    GPR rule analysis and enzyme classification
gene_sequence_mapper
    UniProt sequence retrieval
substrate_mapper
    Metabolite to SMILES mapping
kapp_builder
    Core pipeline functions for kapp calculation
FVA_analysis
    Flux Variability Analysis tools
"""

__version__ = "0.1.0"
