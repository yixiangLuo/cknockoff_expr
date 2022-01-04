#!/bin/bash

N_JOBS=400
EXPR_NAMES=$(./expr_names.sh)

ml R

for EXPR_NAME in $EXPR_NAMES; do
    Rscript ../R/set_expr.R $EXPR_NAME $N_JOBS
done

