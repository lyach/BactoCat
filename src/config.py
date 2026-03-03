"""
BactoCat Configuration Module.

This module provides:
- Static project paths and database file locations
- Pydantic models for pipeline configuration validation
- YAML configuration loading with validation
"""

import os
from pathlib import Path
from typing import Optional, Literal

import yaml
from pydantic import BaseModel, Field, field_validator, model_validator


# =============================================================================
# STATIC PATHS
# =============================================================================

PROJ_ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = PROJ_ROOT / "data"
RAW_DATA_DIR = DATA_DIR / "raw"
INTERIM_DATA_DIR = DATA_DIR / "interim"
PROCESSED_DATA_DIR = DATA_DIR / "processed"
EXTERNAL_DATA_DIR = DATA_DIR / "external_databases"
MODELS_DIR = PROJ_ROOT / "models"
REPORTS_DIR = PROJ_ROOT / "reports"
FIGURES_DIR = REPORTS_DIR / "figures"
RESULTS_DIR = PROJ_ROOT / "results"

# =============================================================================
# DATABASE FILE PATHS
# =============================================================================

BiGG_MAPPING = EXTERNAL_DATA_DIR / "BiGG_mapping.csv"
BIGG_METABOLITES = EXTERNAL_DATA_DIR / "bigg_models_metabolites.txt"
CHEBI_COMPOUNDS = EXTERNAL_DATA_DIR / "CHEBI_compounds.tsv"
CHEBI_INCHI = EXTERNAL_DATA_DIR / "CHEBI_InChI.tsv"
METANETX_COMPOUNDS = EXTERNAL_DATA_DIR / "MetaNetX_compounds.tsv"
METANETX_DEPR = EXTERNAL_DATA_DIR / "MetaNetX_compoundsdepr.tsv"
METANETX_XREF = EXTERNAL_DATA_DIR / "MetaNetX_compoundsxref.tsv"
SEED_COMPOUNDS = EXTERNAL_DATA_DIR / "SEED_compounds.tsv"
SEED_ALIASES = EXTERNAL_DATA_DIR / "Unique_ModelSEED_Compound_Aliases.txt"

# =============================================================================
# ORGANISM CONSTANTS
# =============================================================================

TAXONOMY_IDS = {
    'E coli': 83333,
    'ecoli': 83333,
    'Yeast': 4932,
    'yeast': 4932,
    'S elongatus': 1140,
    'P putida': 160488,
    'B ovatus': 28116,
}

REST_DATASETS = [
    'uniprotkb',
    'uniref50',
    'uniref90',
    'uniref100',
    'tcdb'
]

# =============================================================================
# PIPELINE CONFIGURATION
# =============================================================================


class FluxConfig(BaseModel):
    """Configuration for flux simulation parameters."""
    
    method: Literal["FBA", "pFBA"] = Field(
        default="pFBA",
        description="Flux analysis method: 'FBA' or 'pFBA'"
    )
    carbon_uptake: list[float] = Field(
        description="List of carbon uptake rates to test (mmol/gDW/h)"
    )
    oxygen_uptake: list[float] = Field(
        description="List of oxygen uptake rates to test (mmol/gDW/h)"
    )
    carbon_exchange_rxn: str = Field(
        default="EX_glc__D_e",
        description="Reaction ID for carbon exchange (e.g., glucose uptake)"
    )
    oxygen_exchange_rxn: str = Field(
        default="EX_o2_e",
        description="Reaction ID for oxygen exchange"
    )
    mu_fraction: float = Field(
        default=0.9,
        ge=0.0,
        le=1.0,
        description="Fraction of optimal growth rate for FVA"
    )


class ProteomicsConfig(BaseModel):
    """Configuration for proteomics data parameters."""
    
    p_total: float = Field(
        description="Total protein fraction (g protein / g DCW)"
    )
    paxdb_path: Path = Field(
        description="Path to PaxDB proteomics data file"
    )
    
    @field_validator('paxdb_path', mode='before')
    @classmethod
    def resolve_paxdb_path(cls, v):
        """Convert string to Path and resolve relative paths."""
        if isinstance(v, str):
            return Path(v)
        return v


class ThresholdsConfig(BaseModel):
    """Configuration for kcat filtering thresholds."""
    
    upper_threshold: float = Field(
        default=1.0e6,
        description="Upper limit for kcat_app (s⁻¹)"
    )
    lower_threshold: float = Field(
        default=1.0e-5,
        description="Lower limit for kcat_app (s⁻¹)"
    )


class PipelineConfig(BaseModel):
    """
    Main configuration model for the kapp pipeline.
    
    This model validates YAML configuration files and provides
    type-safe access to all pipeline parameters.
    """
    
    # Model configuration
    organism: str = Field(
        description="Organism name (e.g., 'ecoli', 'yeast')"
    )
    folder_id: str = Field(
        default="run",
        description="ID for output directory organization"
    )
    model_path: Path = Field(
        description="Path to the SBML model file (.xml)"
    )
    solver: str = Field(
        default="cplex",
        description="Solver for optimization (e.g., 'cplex', 'gurobi', 'glpk')"
    )
    
    # Flux simulation
    flux_method: Literal["FBA", "pFBA"] = Field(
        default="pFBA",
        description="Flux analysis method"
    )
    carbon_uptake: Optional[list[float]] = Field(
        default=None,
        description="Carbon uptake rates to test"
    )
    oxygen_uptake: Optional[list[float]] = Field(
        default=None,
        description="Oxygen uptake rates to test"
    )
    carbon_exchange_rxn: str = Field(
        default="EX_glc__D_e",
        description="Carbon exchange reaction ID"
    )
    oxygen_exchange_rxn: str = Field(
        default="EX_o2_e",
        description="Oxygen exchange reaction ID"
    )
    mu_fraction: float = Field(
        default=0.9,
        ge=0.0,
        le=1.0,
        description="Growth rate fraction for FVA"
    )
    medium_df: Optional[Path] = Field(
        default=None,
        description="Path to medium conditions CSV"
    )
    
    # Proteomics
    p_total: float = Field(
        description="Total protein fraction (g protein / g DCW)"
    )
    paxdb_path: Path = Field(
        description="Path to PaxDB file"
    )
    
    # Optional input data
    substrate_df: Optional[Path] = Field(
        default=None,
        description="Path to pre-computed substrate dataframe"
    )
    sequence_df: Optional[Path] = Field(
        default=None,
        description="Path to pre-computed sequence dataframe"
    )
    
    # Thresholds (optional, with defaults)
    upper_threshold: float = Field(
        default=1.0e6,
        description="Upper kcat threshold (s⁻¹)"
    )
    lower_threshold: float = Field(
        default=1.0e-5,
        description="Lower kcat threshold (s⁻¹)"
    )
    calculate_eta: bool = Field(
        default=True,
        description="Whether to compute eta (kapp/kmax) and variance metrics across conditions"
    )
    
    @field_validator('model_path', 'paxdb_path', 'substrate_df', 'sequence_df', 'medium_df', mode='before')
    @classmethod
    def convert_to_path(cls, v):
        """Convert string paths to Path objects."""
        if v is None:
            return None
        if isinstance(v, str):
            return Path(v)
        return v
    
    @model_validator(mode='after')
    def validate_paths_exist(self):
        """Validate that required files exist (when paths are absolute)."""
        # Only validate absolute paths; relative paths are resolved later
        if self.model_path.is_absolute() and not self.model_path.exists():
            raise ValueError(f"Model file not found: {self.model_path}")
        if self.paxdb_path.is_absolute() and not self.paxdb_path.exists():
            raise ValueError(f"PaxDB file not found: {self.paxdb_path}")
        return self
    
    @model_validator(mode='after')
    def validate_flux_inputs(self):
        """Validate that either medium_df, or both carbon/oxygen uptake are provided."""
        if self.medium_df is None:
            if self.carbon_uptake is None or self.oxygen_uptake is None:
                raise ValueError(
                    "Either 'medium_df', or both 'carbon_uptake' and 'oxygen_uptake' must be provided."
                )
        return self
    
    def resolve_paths(self, yaml_dir: Path, project_root: Path) -> "PipelineConfig":
        """
        Resolve relative paths against the project root directory.
        
        Parameters
        ----------
        yaml_dir : Path
            Directory containing the YAML config file
        project_root : Path
            Project root directory
            
        Returns
        -------
        PipelineConfig
            New config with resolved absolute paths
        """
        def resolve(p: Optional[Path]) -> Optional[Path]:
            if p is None:
                return None
            if p.is_absolute():
                return p
            # If path starts with ../, resolve from YAML dir first, then normalize
            if str(p).startswith('../'):
                resolved = (yaml_dir / p).resolve()
                try:
                    resolved.relative_to(project_root)
                    return resolved
                except ValueError:
                    parts = str(p).split('/')
                    up_levels = sum(1 for part in parts if part == '..')
                    remaining = '/'.join(parts[up_levels:])
                    return (project_root / remaining).resolve()
            return (project_root / p).resolve()
        
        return PipelineConfig(
            organism=self.organism,
            folder_id=self.folder_id,
            model_path=resolve(self.model_path),
            solver=self.solver,
            flux_method=self.flux_method,
            carbon_uptake=self.carbon_uptake,
            oxygen_uptake=self.oxygen_uptake,
            carbon_exchange_rxn=self.carbon_exchange_rxn,
            oxygen_exchange_rxn=self.oxygen_exchange_rxn,
            mu_fraction=self.mu_fraction,
            medium_df=resolve(self.medium_df),
            p_total=self.p_total,
            paxdb_path=resolve(self.paxdb_path),
            substrate_df=resolve(self.substrate_df),
            sequence_df=resolve(self.sequence_df),
            upper_threshold=self.upper_threshold,
            lower_threshold=self.lower_threshold,
            calculate_eta=self.calculate_eta,
        )


# =============================================================================
# CONFIGURATION LOADING
# =============================================================================


def load_config(yaml_path: str | Path) -> PipelineConfig:
    """
    Load and validate a pipeline configuration from a YAML file.
    
    Parameters
    ----------
    yaml_path : str or Path
        Path to the YAML configuration file
        
    Returns
    -------
    PipelineConfig
        Validated configuration object with resolved paths
        
    Raises
    ------
    FileNotFoundError
        If the YAML file doesn't exist
    yaml.YAMLError
        If the YAML file is malformed
    pydantic.ValidationError
        If the configuration is invalid
    """
    yaml_path = Path(yaml_path)
    
    if not yaml_path.exists():
        raise FileNotFoundError(f"Configuration file not found: {yaml_path}")
    
    with open(yaml_path, 'r', encoding='utf-8') as f:
        raw_config = yaml.safe_load(f)
    
    # Create config (validates the data)
    config = PipelineConfig(**raw_config)
    
    # Resolve relative paths
    config = config.resolve_paths(yaml_path.parent, PROJ_ROOT)
    
    return config


# =============================================================================
# HELPER FUNCTIONS
# =============================================================================


def ensure_dir_exists(directory: str | Path) -> Path:
    """
    Ensure a directory exists, creating it if necessary.
    
    Parameters
    ----------
    directory : str or Path
        Directory path to create
        
    Returns
    -------
    Path
        The directory path
    """
    directory = Path(directory)
    os.makedirs(directory, exist_ok=True)
    return directory


# =============================================================================
# EXPORTS
# =============================================================================

__all__ = [
    # Paths
    'PROJ_ROOT', 'DATA_DIR', 'RAW_DATA_DIR', 'INTERIM_DATA_DIR',
    'PROCESSED_DATA_DIR', 'EXTERNAL_DATA_DIR', 'MODELS_DIR',
    'REPORTS_DIR', 'FIGURES_DIR', 'RESULTS_DIR',
    # Database files
    'BiGG_MAPPING', 'BIGG_METABOLITES', 'CHEBI_COMPOUNDS', 'CHEBI_INCHI',
    'METANETX_COMPOUNDS', 'METANETX_DEPR', 'METANETX_XREF',
    'SEED_COMPOUNDS', 'SEED_ALIASES',
    # Constants
    'TAXONOMY_IDS', 'REST_DATASETS',
    # Pydantic models
    'PipelineConfig', 'FluxConfig', 'ProteomicsConfig', 'ThresholdsConfig',
    # Functions
    'load_config', 'ensure_dir_exists',
]
