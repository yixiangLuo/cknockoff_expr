#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=31
#SBATCH --nodes=1

#SBATCH --mail-user=yixiangluo@berkeley.edu
#SBATCH --mail-type=ALL

Rscript ../R/do_expr.R main_expr_midScale.R
Rscript ../R/do_expr.R expr_lambdasmax.R
