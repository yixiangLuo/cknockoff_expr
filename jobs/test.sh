#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=7
#SBATCH --nodes=1

#SBATCH --mail-user=yixiangluo@berkeley.edu
#SBATCH --mail-type=ALL

Rscript ../R/do_expr.R test.R
