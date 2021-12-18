#!/bin/bash

Rscript ../R/post_process.R main_expr_largeScale 

Rscript ../R/post_process.R lambdasmax_stat  

Rscript ../R/post_process.R robust_noisyBeta  

Rscript ../R/post_process.R robust_tNoise  

Rscript ../R/post_process.R knockoff_stats  
