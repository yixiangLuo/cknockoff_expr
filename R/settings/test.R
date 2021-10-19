library(here)
library(dbh)
library(knockoff)
library(cknockoff)

library(glmnet)   # lasso
library(KernSmooth)  # local linear regression

source(here("R", "methods.R"))


experiment <- "test"

p <- 1000
n <- 3*p
X_types <- c("IID_Normal", "MCC", "Homo_Block", "Coef_AR", "X_AR")
X_seed <- 2021

pi1 <- 10 / p
posit_types <- c("random", "random", "fix", "random", "random")
beta_permute <- NULL
# beta_permute <- quote(beta[H0] <- runif(sum(H0), min = -0.1/sqrt(n), max = 0.1/sqrt(n)))
noise <- quote(rnorm(n))  # rt(n, df = 5)

target <- 0.5
target_at_alpha <- 0.2
alphas <- c(0.05, 0.1, 0.2, 0.3)

sample_size <- 200
n_cores <- 6

knockoffs <- create.fixed
statistic <- stat.glmnet_coefdiff_lm

get_method_list <- get_multi_method_list
method_names <- c("BH", "BH")

X_titles <- paste0("X: ", X_types)

method_colors <- unname(multi_method_color[method_names])
method_shapes <- unname(multi_method_shape[method_names])











