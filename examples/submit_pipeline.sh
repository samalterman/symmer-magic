#!/bin/bash
#============================================================
# SLURM submission script for the SRE pipeline
#
# Usage:  sbatch submit_pipeline.sh
#============================================================

#SBATCH --job-name=sre-pipeline
#SBATCH --ntasks=32
#SBATCH --mem-per-cpu=4G
#SBATCH --time=02:00:00
#SBATCH --output=sre_pipeline_%j.out
#SBATCH --error=sre_pipeline_%j.err

# Uncomment for GPU jobs:
# #SBATCH --gpus-per-node=4
# #SBATCH --ntasks-per-node=4

#============================================================
# Environment setup
#============================================================

# Load modules (check your cluster)
# module load openmpi/4.1
# module load python/3.11

# Activate conda env
# source activate SRGSP
export PATH="/path/to/your/conda/envs/SRGSP/bin:$PATH"

#============================================================
# Run — edit the Hamiltonian path and options below
#============================================================

HAMILTONIAN="hamiltonian_18_or_less/H2O_STO-3G_SINGLET_JW.json"

echo "Job started: $(date)"
echo "Nodes: $SLURM_JOB_NODELIST"
echo "Tasks: $SLURM_NTASKS"

srun python sre_pipeline.py "$HAMILTONIAN" \
    --solver scipy \
    --k 5 \
    --order 2 \
    --output sre_result.txt

echo "Job finished: $(date)"
