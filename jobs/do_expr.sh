#!/bin/bash

EXPERIMENT=$1
N_JOBS=$2

# set experiment data
ml R
Rscript ../R/set_expr.R $EXPERIMENT $N_JOBS

# submit experiment jobs
sbatch "${EXPERIMENT}_job.sh"
