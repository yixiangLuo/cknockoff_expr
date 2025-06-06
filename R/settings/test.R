library(here)

library(dbh)
library(knockoff)
library(cknockoff)

library(glmnet)   # lasso
library(KernSmooth)  # local linear regression

source(here("R", "utils.R"))
source(here("R", "methods.R"))


experiment <- "test"

p <- 100
n <- 3*p

X_seed <- 2021
pi1 <- 10 / p

X_types <- c("IID_Normal") # "IID_Normal", "MCC", "MCC_Block", "Sparse"
# posit_types <- c("random", "rand_block5")[(X_types %in% c("MCC_Block"))+1]
posit_types <- rep("random", length(X_types))
# random_Xs <- X_types %in% c("IID_Normal")
random_Xs <- rep(T, length(X_types))

# alphas <- c(0.2)
alphas <- c(0.01, 0.05, 0.1, 0.2)
beta_permutes <- NA
# beta_permutes <- c(quote(beta[H0] <- runif(sum(H0), min = -0*mu1, max = 0*mu1)),
#                    quote(beta[H0] <- runif(sum(H0), min = -0.1*mu1, max = 0.1*mu1)),
#                    quote(beta[H0] <- runif(sum(H0), min = -0.2*mu1, max = 0.2*mu1)),
#                    quote(beta[H0] <- runif(sum(H0), min = -0.3*mu1, max = 0.3*mu1)))
noises <- c(quote(rnorm(n)))
# noises <- c(quote(rt(n, df = 9)),
#             quote(rt(n, df = 7)),
#             quote(rt(n, df = 5)),
#             quote(rt(n, df = 3)))
# noises <- c(quote(corr_noise(n, rho = 0)),
#             quote(corr_noise(n, rho = 0.3)),
#             quote(corr_noise(n, rho = 0.6)),
#             quote(corr_noise(n, rho = 0.9)))

fig_x_var <- list(name = "nominal FDR level", value = c(0.01, 0.05, 0.1, 0.2))
# fig_x_var <- list(name = "maximal noise magnitude / mu1", value = c(0, 0.1, 0.2, 0.3))
# fig_x_var <- list(name = "degree of freedom of t-noise", value = c(9, 7, 5, 3))
# fig_x_var <- list(name = "noise equi-correlation rho", value = c(0, 0.3, 0.6, 0.9))
makeup_vectors(alphas = alphas, beta_permutes = beta_permutes, noises = noises)

target <- 0.5
target_at_alpha <- 0.2
calib_method <- "BH" # "BH", "lasso"


sample_size <- 100
n_cores <- 14

knockoffs <- ckn.create.fixed
statistic <- stat.glmnet_coefdiff_tiebreak # ckn.modelX::stat.glmnet_coefdiff_tiebreak

get_method_list <- get_multi_method_list
method_names <- c("BH", "knockoff", "cKnockoff")
# method_names <- c("BH", "dBH", "knockoff", "BonBH", "cKnockoff", "cKnockoff_STAR") #, "mKnockoff"
# method_names <- c("knockoff.MX", "cKnockoff.MX")

X_titles <- paste0("X: ", X_types)

method_colors <- unname(multi_method_color[method_names])
method_shapes <- unname(multi_method_shape[method_names])










