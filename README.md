# 🦠 BactoCat 

A computational framework designed to bridge the gap between *in vitro* and *in vivo* enzyme kinetics. 

By integrating large-scale phenomics, proteomics, and genome-scale metabolic models (GEMs), BactoCat derives apparent catalytic rates ($k_{app}$) that reflect the true catalytic potential of enzymes within bacterial cells.

## Navigation
* [Getting Started](#getting-started)
* [Basic usage](#basic-usage)
* [Repository Structure](#repository-structure)


## Getting Started

1. Clone the repository:
   ```bash
   git clone https://github.com/lyach/BactoCat.git
   ```

2. Navigate to the project directory:
   ```bash
   cd BactoCat
   ```

3. Create an isolated virtual environment:
   
   **Using `uv` (recommended)**
    ```bash
    # Install uv (if first time)
   curl -LsSf https://astral.sh/uv/install.sh | sh # Linux/macOS
   # or
   pip install uv # Windows

   # Create environment
   uv venv --python 3.10 .venv

   # Activate environment
   source .venv/bin/activate # Linux/macOS
   # or 
   .venv\Scripts\activate # Windows

   # Install dependencies
   uv sync
   ```

   **Using `conda`**
    ```bash
    # Create environment
   conda create -n bactocat python=3.10

   # Activate environment
   conda activate bactocat

   # Install dependencies
   pip install -e .
   ```


## Basic Usage

After set-up, you can run the main pipeline  with:
```bash
python scripts/run_kapp_pipeline configs/run_kapp_pipeline/ecoli_medium_aidaying.yaml
```

You can also run a detailed verbose output with `v`:
```bash
python scripts/run_kapp_pipeline -v configs/run_kapp_pipeline/ecoli_medium_aidaying.yaml
```

A full documentation of the pipeline parameters, inputs and outputs can be found [here](docs/kapp_pipeline.md).


## Repository Structure
