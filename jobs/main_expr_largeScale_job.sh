#!/bin/bash

#SBATCH --job-name=expr
#SBATCH --output=../log/main_expr_largeScale/%a.out
#SBATCH --array=1-800
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --time=2-00:00:00

#SBATCH --mail-user=yixiangluo0111@gmail.com
#SBATCH --mail-type=ALL

EXPERIMENT=main_expr_largeScale

# run experiments

ml R

Rscript ../R/experiment.R $EXPERIMENT $SLURM_ARRAY_TASK_ID


