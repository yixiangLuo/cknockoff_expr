#!/usr/bin/env Rscript

library(here)

source(here("R", "experiment.R"))
source(here("R", "methods.R"))
source(here("R", "utils.R"))
source(here("R", "plot.R"))

# source(here("R", "settings", "test.R"))
# source(here("R", "settings", "main_expr_10.R"))
# source(here("R", "settings", "main_expr_30.R"))
# source(here("R", "settings", "lambdasmax_stat.R"))
# source(here("R", "settings", "MVR_kn.R"))
# source(here("R", "settings", "knockoff_stats.R"))
# source(here("R", "settings", "robust_tNoise.R"))
# source(here("R", "settings", "robust_noisyBeta.R"))
# source(here("R", "settings", "cLasso.R"))
source(here("R", "settings", "ModelX.R"))
# source(here("R", "settings", "power_gain.R"))


args <- commandArgs(trailingOnly=TRUE)
if (length(args) > 1) {
  stop("One and only one experiment at a time.")
} else if (length(args) == 1) {
  source(here("R", "settings", paste0(args[1], ".R")))
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
  random_X <- random_Xs[iter_i]
  if(exists("targets")) target <- targets[iter_i]
  if(exists("pi1s")) pi1 <- pi1s[iter_i]
  
  print(paste(X_type, pi1, target))

  if(!random_X){
    X <- gene_X(X_type, n, p, X_seed)$X
    random_X.data <- list(random_X = random_X)
    mu1 <- signal_calib(calib_method, X, random_X.data, pi1, noise = quote(rnorm(n)),
                        posit_type, 1, side = "two", nreps = 20,
                        alpha = target_at_alpha, target = target, n_cores = n_cores)
  } else{
    X.res <- gene_X(X_type, n, p, X_seed)
    X <- NA
    random_X.data <- list(random_X = random_X,
                          n = n, p = p, X_type = X_type, Xcov.true = X.res$Xcov.true,
                          sample_num = 5)
    mu1 <- signal_calib(calib_method, X, random_X.data, pi1, noise = quote(rnorm(n)),
                        posit_type, 1, side = "two", nreps = 20,
                        alpha = target_at_alpha, target = target, n_cores = n_cores)
  }

  problem_setting <- list(n = n, p = p,
                          X = X, X_type = X_type, random_X = random_X,
                          mu1 = mu1, pi1 = pi1, posit_type = posit_type,
                          knockoffs = knockoffs, statistic = statistic,
                          method_names = method_names)

  result <- get_fdp_power(problem_setting, beta_permutes, noises, alphas,
                          fig_x_var, sample_size, n_cores,
                          experiment, X_title)

  save(result, alphas, fig_x_var,
       file = here("data", "temp", paste0(experiment, "-", X_type, ".RData")))

  return(result)
})

names(results) <- X_types

save(results, alphas, fig_x_var, file = here("data", paste0(experiment, ".RData")))
for(X_type in X_types){
  file.remove(here("data", "temp", paste0(experiment, "-", X_type, ".RData")))
}

# method_names <- c("BH", "dBH", "knockoff", "BonBH")
# method_names <- c("BH", "knockoff", "cKnockoff", "cKnockoff_STAR")
# method_colors <- unname(multi_method_color[method_names])
# method_shapes <- unname(multi_method_shape[method_names])
# method_colors <- c("#984ea3", "dodgerblue3", "#333333", "orange1", "red")
# method_shapes <- rep(19, length(method_names))

# X_types <- c("IID_Normal", "MCC", "MCC_Block")
draw_fdp_power_curve(experiment, X_types, sample_size,
                     method_names, method_colors, method_shapes,
                     error_bar = F, direction = F)

# draw_power_gain()

# method_names <- c("BH", "dBH", "knockoff", "cKnockoff", "cKnockoff_STAR")
# method_colors <- unname(multi_method_color[method_names])
# X_types <- c("IID_Normal", "MCC", "MCC_Block")
# 
# draw_fdp_power_dist(experiment, X_types, method_names, method_colors,
#                     type = "FDP", figure = "ECDF")

