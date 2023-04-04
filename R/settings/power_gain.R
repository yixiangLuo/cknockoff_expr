library(here)

library(dbh)
library(knockoff)
library(cknockoff)

library(glmnet)   # lasso
library(KernSmooth)  # local linear regression

source(here("R", "utils.R"))
source(here("R", "methods.R"))


experiment <- "power_gain"

p <- 500
n <- 3*p

X_seed <- 2021
pi1s <- c(5, 10, 20, 30) / p

X_types <- c("IID_Normal", "MCC", "MCC_Block")
posit_types <- rep("random", length(X_types))
random_Xs <- rep(T, length(X_types))

alphas <- c(0.05)
beta_permutes <- NA
noises <- c(quote(rnorm(n)))

fig_x_var <- list(name = "nominal FDR level", value = c(0.05))
makeup_vectors(alphas = alphas, beta_permutes = beta_permutes, noises = noises)

targets <- c(0.5, 0.8, 0.95)
target_at_alpha <- 0.2
calib_method <- "BH" # "BH", "lasso"

x_len <- length(pi1s)
y_len <- length(targets)
set_len <- length(X_types)
X_types <- c(sapply(X_types, function(X_type){
  tmp <- c()
  for(i in 1:y_len){
    for (j in 1:x_len){
      tmp = c(tmp, paste(X_type, targets[i], pi1s[j]*p, sep = "_D_"))
    }
  }
  return(tmp)
}))
posit_types <- rep("random", length(X_types))
random_Xs <- rep(T, length(X_types))
pi1s <- rep(rep(pi1s, y_len), set_len)
targets <- rep(rep(targets, each = x_len), set_len)

sample_size <- 100
n_cores <- 14

knockoffs <- ckn.create.fixed
statistic <- stat.glmnet_coefdiff_tiebreak

get_method_list <- get_multi_method_list
method_names <- c("knockoff", "cKnockoff")

X_titles <- paste0("X: ", X_types)

method_colors <- unname(multi_method_color[method_names])
method_shapes <- unname(multi_method_shape[method_names])










