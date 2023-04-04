library(here)
library(glmnet) 

library(tidyverse)
library(ggpubr) 
library(latex2exp)

library(foreach)
library(doParallel)


source(here("R", "utils.R"))
source(here("R", "methods.R"))


experiment <- "ROC_compare"

p <- 1000
n <- 3*p
X_type <- "MCC_Block"
X_seed <- 2021

pi1 <- 10 / p
posit_type <- "random"
noise <- quote(rnorm(n))

target <- 0.5
target_at_alpha <- 0.2

sample_size <- 100

n_cores <- 14


X <- gene_X(X_type, n, p, X_seed)$X

random_X.data <- list(random_X = F)
mu1 <- BH_lm_calib(X, random_X.data, pi1, noise, posit_type, 1,
                   side = "two", nreps = 200,
                   alpha = target_at_alpha, target = target, n_cores = n_cores)
beta <- genmu(p, pi1, mu1, posit_type, 1)

H0 <- beta == 0


registerDoParallel(n_cores)


results <- foreach(iter = 1:sample_size) %dopar% {
    
    set.seed(iter)
    
    y <- X %*% beta + eval(noise)
    
    OLS_fit <- lm(y ~ X + 0)
    OLS_beta <- abs(OLS_fit$coefficients)
    
    sigma_hat <- sqrt(sum(OLS_fit$residuals^2) / (n-p))
    lambda <-  2 * sigma_hat / n
    LASSO_fit <- glmnet::glmnet(X, y, lambda = lambda,
                                intercept = F, standardize = F,
                                standardize.response = F, family = "gaussian")
    LASSO_beta <- abs(LASSO_fit$beta)
    y_res <- as.matrix(y - (X %*% LASSO_fit$beta))
    LM <- c(abs(t(y_res) %*% X))
    LM_add_on <- rep(lambda * n, p)
    not_in_model <- which(LASSO_beta < 1e-10)
    LM_add_on[not_in_model] <- LM[not_in_model]
    LASSO_beta <- LASSO_beta + LM_add_on
    
    OLS_ROC <- sapply(sort(OLS_beta), function(thr){
        selected <- which(OLS_beta >= thr)
        FDP <- length(intersect(selected, which(H0))) / max(length(selected), 1)
        FPP <- length(intersect(selected, which(H0))) / ((1-pi1) * p)
        TPP <- length(intersect(selected, which(!H0))) / (pi1 * p)
        return(c(FPP, FDP, TPP, length(selected)))
    })
    OLS_ROC <- data.frame(FPP = OLS_ROC[1, ], FDP = OLS_ROC[2, ],
                          TPP = OLS_ROC[3, ], n_sel = OLS_ROC[4, ],
                          stat = "OLS") %>%
        mutate(FPP = FPP + 1e-10 * TPP, FDP = FDP + 1e-10 * TPP)
    
    LASSO_ROC <- sapply(sort(LASSO_beta), function(thr){
        selected <- which(LASSO_beta >= thr)
        FDP <- length(intersect(selected, which(H0))) / max(length(selected), 1)
        FPP <- length(intersect(selected, which(H0))) / ((1-pi1) * p)
        TPP <- length(intersect(selected, which(!H0))) / (pi1 * p)
        return(c(FPP, FDP, TPP, length(selected)))
    })
    LASSO_ROC <- data.frame(FPP = LASSO_ROC[1, ], FDP = LASSO_ROC[2, ],
                            TPP = LASSO_ROC[3, ], n_sel = LASSO_ROC[4, ],
                            stat = "LASSO") %>%
        mutate(FPP = FPP + 1e-10 * TPP, FDP = FDP + 1e-10 * TPP)
    
    return(list(OLS_ROC = OLS_ROC, LASSO_ROC = LASSO_ROC))
}

grid_num <- 5000
x_points <- seq(from = 0, to = 0.5, length.out = grid_num)
x_var <- "FPP" # FDP

TPP_OLS <- sapply(results, function(result){
    approx(x = result$OLS_ROC[[x_var]], 
           y = result$OLS_ROC$TPP, 
           xout = x_points,
           method="linear")$y
})
TPP_OLS <- rowMeans(TPP_OLS, na.rm = T)

TPP_LASSO <- sapply(results, function(result){
    approx(x = result$LASSO_ROC[[x_var]], 
           y = result$LASSO_ROC$TPP, 
           xout = x_points,
           method="linear")$y
})
TPP_LASSO <- rowMeans(TPP_LASSO, na.rm = T)

plot_data <- data.frame(var = rep(x_points, 2),
                        TPP = c(TPP_OLS, TPP_LASSO),
                        stat = rep(c("OLS", "LASSO"), each = grid_num))

method_colors <- c("dodgerblue3", "red")
# method_names <- unname(TeX(c("$|\\hat{\\beta}_{OLS}|$", "$|\\hat{\\beta}_{LASSO}|$")))
method_names <- c("OLS", "LASSO")

ggplot(plot_data) +
    geom_line(aes(x = var, y = TPP, color = stat)) +
    scale_color_manual(values = method_colors, breaks = c("OLS", "LASSO"),
                       labels = method_names) +
    scale_x_continuous(limits = c(0, 0.5), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    theme_bw() +
    theme(aspect.ratio = 1,
          # panel.grid = element_blank(),
          panel.border = element_blank(),
          strip.text = element_text(size = 13),
          axis.title = element_text(size = 11),
          axis.text = element_text(size = 8),
          legend.position = "right",
          legend.title=element_text(size=9),
          legend.text=element_text(size=9)) +
    labs(color = "variable ordering", x = x_var)

# figure <- ggarrange(plot1, plot2,
#                     # labels = c("A", "B", "C"),
#                     ncol = 2, nrow = 1,
#                     common.legend = T, legend = "right")

ggsave(filename = here("figs", paste0("ROC_compare.pdf")),
       width = 4, height = 3)
# ggsave(filename = here("figs", paste0("ROC_compare.pdf")),
#        figure, width = 6.5, height = 2.5)


