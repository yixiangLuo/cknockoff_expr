library(here)

library(dbh)
library(knockoff)
library(cknockoff)

library(glmnet)   # lasso
library(KernSmooth)  # local linear regression

source(here("R", "utils.R"))
source(here("R", "methods.R"))

# hard-coded settings
# 1. utils.R gene_X: model_X <- T # set False if experiment with fixed-X
# 2. utils.R MXkn_calib: family = "gaussain" / "binomial"
# 3. experiment.R: y generated as 0/1 or gaussian
# 4. do_expr.R/set_expr.R: signal_calib nreps = 20 / 200


experiment <- "ModelX"

p <- 100
n <- 300

X_seed <- 2021
pi1 <- 10 / p

X_types <- c("IID_Normal", "MCC", "MCC_Block")
posit_types <- rep("random", length(X_types))
random_Xs <- rep(T, length(X_types))

alphas <- c(0.01, 0.05, 0.1, 0.2)
beta_permutes <- NA
noises <- c(quote(rnorm(n)))

fig_x_var <- list(name = "nominal FDR level", value = c(0.01, 0.05, 0.1, 0.2))
makeup_vectors(alphas = alphas, beta_permutes = beta_permutes, noises = noises)

targets <- c(0.6, 0.1, 0.7)
target_at_alpha <- 0.2
calib_method <- "MXkn" # "BH", "lasso", "MXkn"


sample_size <- 100
n_cores <- 14

knockoffs <- ckn.create.fixed
# statistic <- ckn.modelX::stat.glmnet_coefdiff_tiebreak
# statistic <- ckn.modelX::stat.glmnet_coefdiff_cv
statistic <- function(X, X_k, y){
  ckn.modelX::stat.glmnet_coefdiff_tiebreak(X, X_k, y, sigma_tilde = 1,
                                            family = "binomial")
}

get_method_list <- get_multi_method_list
method_names <- c("knockoff.MX", "cKnockoff.MX")

X_titles <- paste0("X: ", X_types)

method_colors <- unname(multi_method_color[method_names])
method_shapes <- unname(multi_method_shape[method_names])










