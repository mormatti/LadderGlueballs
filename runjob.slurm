#!/bin/bash

#SBATCH --job-name=glueballs
#SBATCH --output=job_out/glueballs.out.job-%j
#SBATCH --error=job_out/glueballs.err.job-%j
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=25
#SBATCH --partition=long
#SBATCH --time=30-00:00:00

cd /mnt/mmorgavi/LadderGlueballs

# Usa il progetto locale
export JULIA_PROJECT=.

# Dillo a Julia e a BLAS quanti core hai DAVVERO
export JULIA_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK   # nel dubbio, non fa male

julia --color=yes --project=. 4_creator.jl | tee job_out/live.log