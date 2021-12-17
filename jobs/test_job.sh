#!/bin/bash

#SBATCH --job-name=expr
#SBATCH --output=../log/test/%a.out
#SBATCH --array=1-5
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=2-00:00:00

EXPERIMENT=test

# run experiments
Rscript ../R/experiment.R $EXPERIMENT $SLURM_ARRAY_TASK_ID


