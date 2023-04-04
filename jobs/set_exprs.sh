#!/bin/bash

#SBATCH --job-name=set_expr
#SBATCH --output=../log/%a.out
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --time=1-00:00

#SBATCH --mail-user=email@address
#SBATCH --mail-type=ALL

N_JOBS=400
EXPR_NAMES=$(./expr_names.sh)

ml R

for EXPR_NAME in $EXPR_NAMES; do
    Rscript ../R/cluster/set_expr.R $EXPR_NAME $N_JOBS
done

