#!/bin/bash
#SBATCH -p batch
#SBATCH -t 0-12:00:00
#SBATCH --constraint=[ivb|hsw]
#SBATCH -n 1
#SBATCH --mem=3000
#SBATCH --array=1-30
#SBATCH -o out/run-%A-%a.out
module load r/3.6.1-python3
srun Rscript --vanilla run.R $SLURM_ARRAY_TASK_ID
