#!/bin/bash
#SBATCH --job-name=RnaseqSnake
#SBATCH --output=miworkflow_%j.log
#SBATCH --error=miworkflow_%j.err
#SBATCH -p cicese
#SBATCH --ntasks-per-node=24
#SBATCH --mem=100GB
#SBATCH -t 6-00:00:00  # Cambia esto según lo que necesites

# Cargar Conda/Mamba
source /LUSTRE/apps/Miniforge/2024/miniforge3/etc/profile.d/conda.sh
source /LUSTRE/apps/Miniforge/2024/miniforge3/etc/profile.d/mamba.sh
CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh
conda activate miworkflow_env  # cambia al entorno que uses

cd $SLURM_SUBMIT_DIR

snakemake -R $(snakemake --list-params-changes) \
 -d ./ \
 --snakefile snakefile \
 --workflow-profile ./profiles/slurm \
 --cores 24 \
 --use-conda
