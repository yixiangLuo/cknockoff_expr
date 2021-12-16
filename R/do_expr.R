#!/usr/bin/env Rscript

library(here)

source(here("R", "experiment.R"))
source(here("R", "methods.R"))
source(here("R", "utils.R"))
source(here("R", "plot.R"))

# source(here("R", "settings", "test.R"))
# source(here("R", "settings", "main_expr_smallScale.R"))
# source(here("R", "settings", "main_expr_midScale.R"))
# source(here("R", "settings", "main_expr_largeScale.R"))
# source(here("R", "settings", "knockoff_stats.R"))

args <- commandArgs(trailingOnly=TRUE)
if (length(args) > 1) {
  stop("One and only one experiment at a time.")
} else if (length(args) == 1) {
  source(here("R", "settings", args[1]))
}

file_name <- here("data", "temp", paste0("progress-", experiment, ".txt"))
if (file.exists(file_name)) {
  file.remove(file_name)
  invisible()
}


results <- lapply(1:length(X_types), function(iter_i){
  X_type <- X_types[iter_i]
  posit_type <- posit_types[iter_i]
  X_title <- X_titles[iter_i]
  
  X <- gene_X(X_type, n, p, X_seed)
  
  mu1 <- BH_lm_calib(X, pi1, noise = quote(rnorm(n)),
                     posit_type, 1, side = "two", nreps = 200,
                     alpha = target_at_alpha, target = target, n_cores = n_cores)
  beta <- genmu(p, pi1, mu1, posit_type, 1)
  
  H0 <- beta == 0
  
  method_list <- get_method_list(X, knockoffs, statistic)[method_names]
  
  result <- get_fdp_power(X, beta, H0, mu1, beta_permutes, noises, alphas,
                          fig_x_var, method_list,
                          sample_size, n_cores,
                          experiment, X_title)
  
  save(result, alphas, fig_x_var,
       file = here("data", "temp", paste0(experiment, "-", X_type, ".Rdata")))
  
  return(result)
})

names(results) <- X_types

save(results, alphas, fig_x_var, file = here("data", paste0(experiment, ".Rdata")))
for(X_type in X_types){
  file.remove(here("data", "temp", paste0(experiment, "-", X_type, ".Rdata")))
}

# method_names <- c("BH", "dBH", "knockoff", "BonBH")
# method_colors <- unname(multi_method_color[method_names])
# method_shapes <- unname(multi_method_shape[method_names])
# method_shapes <- rep(19, 4)

# X_types <- c("IID_Normal", "MCC", "Homo_Block")
draw_fdp_power_curve(experiment, X_types, sample_size,
                     method_names, method_colors, method_shapes,
                     error_bar = F, direction = F)

