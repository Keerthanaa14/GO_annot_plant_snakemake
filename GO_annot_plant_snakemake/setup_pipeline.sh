#!/bin/bash
#SBATCH --account=project_200xxx   # your Puhti project ID
#SBATCH --time=03:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --job-name=setup_annotation
#SBATCH --output=logs/setup_%j.log

# Exit on error
set -euo pipefail

# Load required modules on Puhti
module load biokit
module load blast
module load diamond
module load python/3.9

# === 1. Create project structure ===
mkdir -p ~/plant_annotation/{data,results,envs,logs,scripts}
cd ~/plant_annotation

# === 2. Extract plant sequences from nr ===
echo "Extracting plant protein sequences from nr database..."
cd data

# List of taxon IDs for plants (Viridiplantae and subgroups)
# Viridiplantae (33090) includes all green plants
PLANT_TAXIDS=33090

blastdbcmd \
  -taxids ${PLANT_TAXIDS} \
  -db nr \
  -dbtype prot \
  -out nr_plants.fasta \
  -target_only

echo "✅ Plant FASTA extracted: data/nr_plants.fasta"

# === 3. Convert to DIAMOND database ===
echo "Converting to DIAMOND database..."
diamond makedb \
  --in nr_plants.fasta \
  --db nr_plants.dmnd \
  --threads 8

echo "✅ DIAMOND database created: data/nr_plants.dmnd"

# === 4. Clone EggNOG-mapper ===
cd ~/plant_annotation
echo "Cloning EggNOG-mapper repository..."
git clone https://github.com/eggnogdb/eggnog-mapper.git
cd eggnog-mapper

# === 5. Download EggNOG database files ===
echo "Downloading EggNOG database..."
mkdir -p data
cd data

BASE_URL="http://eggnog5.embl.de/download/emapperdb-5.0.2"

wget -c ${BASE_URL}/eggnog.db.gz
wget -c ${BASE_URL}/eggnog.taxa.tar.gz
wget -c ${BASE_URL}/eggnog_proteins.dmnd.gz
wget -c ${BASE_URL}/eggnog_proteins_hmm.tar.gz

# Unpack all files
gunzip eggnog.db.gz
gunzip eggnog_proteins.dmnd.gz
tar -xzf eggnog.taxa.tar.gz
tar -xzf eggnog_proteins_hmm.tar.gz

echo "✅ EggNOG database prepared under: eggnog-mapper/data"

# === 6. Create a Python virtual environment (no conda) ===
cd ~/plant_annotation
python3 -m venv envs/eggnog
source envs/eggnog/bin/activate

pip install --upgrade pip
pip install -e ./eggnog-mapper

echo "✅ EggNOG-mapper environment ready and activated."

# === 7. Verify installation ===
echo "Verifying EggNOG-mapper setup..."
python -m emapper --help

echo "🎉 Setup complete! You can now use your custom plant DIAMOND database with EggNOG-mapper."