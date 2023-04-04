#!/bin/bash

#SBATCH --job-name=power_gain
#SBATCH --output=../log/power_gain/%a.out
#SBATCH --array=1-14
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --time=2-00:00:00

#SBATCH --mail-user=email@address
#SBATCH --mail-type=ALL

EXPERIMENT=power_gain

# run experiments

ml R

Rscript ../R/cluster/experiment.R $EXPERIMENT $SLURM_ARRAY_TASK_ID


