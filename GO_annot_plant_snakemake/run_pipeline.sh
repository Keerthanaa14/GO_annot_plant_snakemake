#!/bin/bash
#SBATCH --job-name=plant_annot
#SBATCH --time=12:00:00
#SBATCH --mem=512G
#SBATCH --cpus-per-task=32
#SBATCH --output=logs/snakemake_%j.out

set -euo pipefail
mkdir -p logs

echo "Starting workflow at $(date '+%Y-%m-%d_%H-%M-%S')"

# Use SLURM_CPUS_PER_TASK if available, otherwise default to 1
CPUS=${SLURM_CPUS_PER_TASK:-32}

# Load required modules
module load python-data
module load snakemake
module load biokit
module load diamond
module load interproscan

# Run Snakemake workflow
snakemake \
    --snakefile ~/GO_annot_plant_snakemake/Snakefile \
    --configfile ~/GO_annot_plant_snakemake/config.yaml \
    --cores "$CPUS" \
    --rerun-incomplete \
    --printshellcmds \
    --latency-wait 60 \
    --keep-going

echo "Finished workflow at $(date '+%Y-%m-%d_%H-%M-%S')"