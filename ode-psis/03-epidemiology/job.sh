#!/bin/bash
#SBATCH -p batch
#SBATCH -t 0-04:00:00
#SBATCH --array=1
#SBATCH --mem=3000
#SBATCH -o out/out-%a.out
module load r/4.0.3-python3
srun Rscript --vanilla sample_rk45.R $SLURM_ARRAY_TASK_ID
