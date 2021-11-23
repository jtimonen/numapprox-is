#!/bin/bash
#SBATCH -p batch
#SBATCH -t 2-08:00:00
#SBATCH --array=1-60
#SBATCH --mem=3000
#SBATCH -o out/job-%a.out
module load r/4.0.3-python3
srun Rscript --vanilla main.R $SLURM_ARRAY_TASK_ID
