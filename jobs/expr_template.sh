#!/bin/bash

#SBATCH --job-name=expr
#SBATCH --output=../log/expr/%a.out
#SBATCH --array=1-N_JOBS
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --time=2-00:00:00

#SBATCH --mail-user=yixiangluo0111@gmail.com
#SBATCH --mail-type=ALL

EXPERIMENT=expr

# run experiments
Rscript ../R/experiment.R $EXPERIMENT $SLURM_ARRAY_TASK_ID


