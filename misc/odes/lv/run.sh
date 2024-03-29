#!/bin/bash
#SBATCH -p batch
#SBATCH -t 0-12:00:00
#SBATCH --constraint=[ivb|hsw]
#SBATCH -n 1
#SBATCH --mem=5000
#SBATCH --array=1-10
#SBATCH -o out/rk4/run-%A-%a.out
module load r/3.6.1-python3
srun Rscript --vanilla run_rk4.R $SLURM_ARRAY_TASK_ID
