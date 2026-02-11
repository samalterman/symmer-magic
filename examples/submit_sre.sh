#!/bin/bash
#============================================================
# SLURM submission script for SRE computation with MPI
# Usage:  sbatch submit_sre.sh
#============================================================

#--- Job name ---
#SBATCH --job-name=sre-calc

#--- Number of MPI tasks (= number of CPU cores doing work) ---
#--- Increase this for larger qubit counts ---
#SBATCH --ntasks=32

#--- How many nodes to spread across (optional, SLURM can auto-decide) ---
# #SBATCH --nodes=2

#--- Memory per CPU core ---
#SBATCH --mem-per-cpu=4G

#--- Wall clock time limit (HH:MM:SS) ---
#--- For 14-16 qubits with 32 ranks, ~1 hour is usually enough ---
#SBATCH --time=02:00:00

#--- Output and error log files (%j = job ID) ---
#SBATCH --output=sre_%j.out
#SBATCH --error=sre_%j.err

#--- Partition / queue name (check cluster docs) ---
# #SBATCH --partition=compute
# #SBATCH --partition=gpu        # uncomment for GPU jobs

#--- GPU options (uncomment for MPI+GPU jobs) ---
# #SBATCH --gpus-per-node=4
# #SBATCH --ntasks-per-node=4    # typically 1 MPI rank per GPU

# Load modules 
#*****************

# Activate your conda environment
#******************
#============================================================
# Run
#============================================================

echo "Job started: $(date)"
echo "Nodes: $SLURM_JOB_NODELIST"
echo "Tasks: $SLURM_NTASKS"

# srun is preferred over mpirun on SLURM clusters
srun python run_sre_mpi.py

echo "Job finished: $(date)"
