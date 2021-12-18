#!/bin/bash

N_JOBS=800

sh do_expr.sh main_expr_largeScale $N_JOBS

sh do_expr.sh lambdasmax_stat $N_JOBS

sh do_expr.sh robust_noisyBeta $N_JOBS

sh do_expr.sh robust_tNoise $N_JOBS

sh do_expr.sh knockoff_stats $N_JOBS
