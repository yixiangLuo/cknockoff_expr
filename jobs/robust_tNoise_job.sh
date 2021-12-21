#!/bin/bash

#SBATCH --job-name=expr
#SBATCH --output=../log/robust_tNoise/%a.out
#SBATCH --array=1-400
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --time=2-00:00:00

#SBATCH --mail-user=yixiangluo0111@gmail.com
#SBATCH --mail-type=ALL

EXPERIMENT=robust_tNoise

# run experiments

ml R

Rscript ../R/experiment.R $EXPERIMENT $SLURM_ARRAY_TASK_ID


