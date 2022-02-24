#!/bin/bash

#SBATCH --job-name=main_expr_30
#SBATCH --output=../log/main_expr_30/%a.out
#SBATCH --array=1-400
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --time=2-00:00:00

#SBATCH --mail-user=yixiangluo@berkeley.edu
#SBATCH --mail-type=ALL

EXPERIMENT=main_expr_30

# run experiments

ml R

Rscript ../R/experiment.R $EXPERIMENT $SLURM_ARRAY_TASK_ID


