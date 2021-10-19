
library(knockoff)
library(cknockoff)

library(glmnet)   # lasso

source(here("R", "methods.R"))


experiment <- "knockoff_stats"

p <- 300
n <- 3*(p+1)
X_types <- c("IID_Normal", "MCC", "Homo_Block", "Coef_AR", "X_AR")
X_seed <- 2021

pi1 <- 10 / p
posit_types <- c("random", "random", "fix", "random")
beta_permute <- NULL
noise <- quote(rnorm(n))

target <- 0.5
target_at_alpha <- 0.2
alphas <- c(0.05, 0.1, 0.2, 0.3)

sample_size <- 300
n_cores <- 6

knockoffs <- create.fixed
statistic <- stat.glmnet_coefdiff_lm

get_method_list <- get_kn_method_list
method_names <- c("kn_D_lambdasmax", "kn_D_lambdasmax_lm", "kn_D_coefdiff_lm")

X_titles <- paste0("X: ", X_types)

method_colors <- unname(kn_method_color[method_names])
method_shapes <- unname(kn_method_shape[method_names])











