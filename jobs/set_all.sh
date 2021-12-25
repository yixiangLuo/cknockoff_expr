#!/bin/bash

N_JOBS=400

sh set_expr.sh main_expr_largeScale $N_JOBS

sh set_expr.sh lambdasmax_stat $N_JOBS

sh set_expr.sh robust_noisyBeta $N_JOBS

sh set_expr.sh robust_tNoise $N_JOBS

sh set_expr.sh knockoff_stats $N_JOBS
