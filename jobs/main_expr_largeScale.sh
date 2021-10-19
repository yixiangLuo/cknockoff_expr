#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20

Rscript ../R/do_expr.R main_expr_largeScale.R
