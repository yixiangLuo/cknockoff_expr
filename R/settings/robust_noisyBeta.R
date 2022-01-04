library(here)

library(dbh)
library(knockoff)
library(cknockoff)

library(glmnet)   # lasso
library(KernSmooth)  # local linear regression

source(here("R", "utils.R"))
source(here("R", "methods.R"))


experiment <- "robust_noisyBeta"

p <- 1000
n <- 3*p

X_seed <- 2021
pi1 <- 10 / p

X_types <- c("IID_Normal", "MCC", "MCC_Block")
posit_types <- rep("random", length(X_types))

alphas <- c(0.2)
beta_permutes <- c(quote(beta[H0] <- runif(sum(H0), min = -0*mu1, max = 0*mu1)),
                   quote(beta[H0] <- runif(sum(H0), min = -0.1*mu1, max = 0.1*mu1)),
                   quote(beta[H0] <- runif(sum(H0), min = -0.2*mu1, max = 0.2*mu1)),
                   quote(beta[H0] <- runif(sum(H0), min = -0.3*mu1, max = 0.3*mu1)))
noises <- c(quote(rnorm(n)))

fig_x_var <- list(name = TeX("maximal noise magnitude / $\\beta^*$"), value = c(0, 0.1, 0.2, 0.3))
makeup_vectors(alphas = alphas, beta_permutes = beta_permutes, noises = noises)

target <- 0.5
target_at_alpha <- 0.2

sample_size <- 400
n_cores <- 14

knockoffs <- create.fixed
statistic <- stat.glmnet_coefdiff_lm

get_method_list <- get_multi_method_list
method_names <- c("BH", "dBH", "knockoff", "BonBH", "cKnockoff", "cKnockoff_STAR")

X_titles <- paste0("X: ", X_types)

method_colors <- unname(multi_method_color[method_names])
method_shapes <- unname(multi_method_shape[method_names])










