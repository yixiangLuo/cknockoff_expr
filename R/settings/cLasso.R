library(here)

library(dbh)
library(knockoff)
library(cknockoff)

library(glmnet)   # lasso
library(KernSmooth)  # local linear regression

source(here("R", "utils.R"))
source(here("R", "methods.R"))
source(here("R", "cLasso", "cLasso.R"))


experiment <- "cLasso"

p <- 100
n <- 3*p

X_seed <- 2021
pi1 <- 10 / p

X_types <- c("IID_Normal", "MCC_Block", "Coef_AR", "X_AR") # "IID_Normal", "MCC", "MCC_Block", "Sparse"
posit_types <- rep("random", length(X_types))
random_Xs <- X_types %in% c("IID_Normal")

# alphas <- c(0.2)
alphas <- c(0.01, 0.05, 0.1, 0.2)
beta_permutes <- NA
noises <- c(quote(rnorm(n)))

fig_x_var <- list(name = "nominal FDR level", value = c(0.01, 0.05, 0.1, 0.2))
makeup_vectors(alphas = alphas, beta_permutes = beta_permutes, noises = noises)

target <- 0.5
target_at_alpha <- 0.2
calib_method <- "BH"

sample_size <- 400
n_cores <- 14

knockoffs <- ckn.create.fixed
statistic <- stat.glmnet_coefdiff_tiebreak

get_method_list <- get_multi_method_list
method_names <- c("BH", "dBH", "knockoff", "cKnockoff", "cLasso") 

X_titles <- paste0("X: ", X_types)

method_colors <- unname(multi_method_color[method_names])
method_shapes <- unname(multi_method_shape[method_names])










