#!/bin/bash
#SBATCH --job-name=fission_sim
#SBATCH --array=1-400
#SBATCH --time=02:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/out_%A_%a.txt
#SBATCH --error=logs/err_%A_%a.txt

module load R

mkdir -p logs
mkdir -p results

Rscript run_chunk.R ${SLURM_ARRAY_TASK_ID}
