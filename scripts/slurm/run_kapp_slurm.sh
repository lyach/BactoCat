#!/bin/bash
#SBATCH --job-name=kapp_pipeline
#SBATCH --output=logs/kapp_%j.out
#SBATCH --error=logs/kapp_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yuxii.zhang@mail.utoronto.ca

# --- Start info ---
echo "=============================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo "Start time: $(date)"
echo "=============================================="

# --- Project Setup ---
# Set the root of the BactoCat pipeline
PROJECT_DIR="/home/yuxiz/BactoCat" 
cd "$PROJECT_DIR"

# Ensure log directory exists so Slurm doesn't fail
mkdir -p logs

# --- Environment ---
# Activate the virtual environment 
source "$PROJECT_DIR/.venv/bin/activate"

# --- Run Command ---
python -m scripts.run_kapp_pipeline \
    configs/run_kapp_pipeline/ecoli_medium_test.yaml 

# --- End info ---
echo "=============================================="
echo "End time: $(date)"
echo "=============================================="