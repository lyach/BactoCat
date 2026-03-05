"""
BactoCat ETA analysis Runner.

This script runs the ETA analysis for .

Usage:
    run-eta-analysis config.yaml
    python -m scripts.run_eta_analysis config.yaml
"""

import argparse
from pathlib import Path
import pandas as pd
import yaml
from src.config import EtaInVitroConfig, PROJ_ROOT
from src.utils import get_eta_in_vitro

def main():
    parser = argparse.ArgumentParser(description="Run Eta In Vitro Analysis")
    parser.add_argument("config", type=str, help="Path to YAML config")
    args = parser.parse_args()

    # 1. Load raw YAML
    with open(args.config, 'r') as f:
        raw_dict = yaml.safe_load(f)

    # 2. Validate with Pydantic (using the model we added to config.py)
    # This automatically handles the Path conversion
    cfg = EtaInVitroConfig(**raw_dict)

    # 3. Resolve paths relative to Project Root (Standardizing with your setup)
    full_kmax_path = (PROJ_ROOT / cfg.kmax_path).resolve()
    full_in_vitro_path = (PROJ_ROOT / cfg.in_vitro_kcat_path).resolve()

    # 4. Execute
    kmax_df = pd.read_csv(full_kmax_path)
    
    get_eta_in_vitro(
        in_vitro_kcat_path=str(full_in_vitro_path),
        kmax_results=kmax_df,
        kmax_path=str(full_kmax_path),
        dataset_name=cfg.dataset_name
    )

if __name__ == "__main__":
    main()