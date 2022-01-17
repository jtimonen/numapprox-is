#!/bin/bash
#SBATCH -p batch
#SBATCH -t 0-02:00:00
#SBATCH --array=1-60
#SBATCH --mem=3000
#SBATCH -o out/run-%a.out
module load r/4.0.3-python3
srun Rscript --vanilla main_psis.R $SLURM_ARRAY_TASK_ID