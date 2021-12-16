library(here)

library(dbh)
library(knockoff)
library(cknockoff)

library(glmnet)   # lasso
library(KernSmooth)  # local linear regression

source(here("R", "methods.R"))


experiment <- "main_expr_midScale"

p <- 300
n <- 3*p
X_types <- c("IID_Normal", "MCC", "Homo_Block", "MCC_Block")
# X_types <- c("IID_Normal", "MCC", "Homo_Block", "Coef_AR", "X_AR")
X_seed <- 2021

pi1 <- 10 / p
posit_types <- c("random", "random", "fix", "random")
# posit_types <- c("random", "random", "fix", "random", "random")
beta_permute <- NULL
noise <- quote(rnorm(n))

target <- 0.5
target_at_alpha <- 0.2
alphas <- c(0.05, 0.1, 0.2, 0.3)

# sample_size <- 1000
sample_size <- 100
n_cores <- 7

knockoffs <- create.fixed
statistic <- stat.glmnet_coefdiff_lm

get_method_list <- get_multi_method_list
# method_names <- c("BH", "BonBH")
method_names <- c("BH", "dBH", "knockoff", "BonBH", "cKnockoff", "cKnockoff_PLUS")

X_titles <- paste0("X: ", X_types)

method_colors <- unname(multi_method_color[method_names])
method_shapes <- unname(multi_method_shape[method_names])





# for(X_type in X_types){
#   BBH <- results_BBH[[X_type]]$FDR_Power %>% filter(method_name == "BonBH") %>%
#     mutate(methods = factor(4, levels = 1:6))
#   FDR_Power <- rbind(results_org[[X_type]]$FDR_Power %>% filter(method_name != "BonBH"),
#                      BBH)
#   results[[X_type]]$FDR_Power <- FDR_Power
# }
# save(results, file = here("data", paste0("main_expr_midScale_adj", ".Rdata")))




