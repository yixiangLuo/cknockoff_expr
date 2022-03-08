#!/bin/bash

EXPR_NAMES=$(./expr_names.sh)

ml R

for EXPR_NAME in $EXPR_NAMES; do
    Rscript ../R/post_process.R $EXPR_NAME 
done
