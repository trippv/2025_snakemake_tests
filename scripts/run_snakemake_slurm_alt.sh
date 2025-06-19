#!/bin/bash
#SBATCH --job-name=RnaseqSnake
#SBATCH --output=rnaseq_snake_%j.log
#SBATCH --error=rnaseq_snake_%j.err
#SBATCH -p cicese
#SBATCH --ntasks-per-node=24
#SBATCH --mem=100GB
#SBATCH -t 6-00:00:00  # Cambia esto según lo que necesites

# Cargar Conda/Mamba
source /LUSTRE/apps/Miniforge/2024/miniforge3/etc/profile.d/conda.sh


conda activate snakemake  # cambia al entorno que uses

cd $SLURM_SUBMIT_DIR


# Ejecutar Snakemake en modo local con 24 núcleos
snakemake --cores 24 --use-conda

