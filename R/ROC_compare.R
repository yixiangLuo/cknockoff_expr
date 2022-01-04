library(here)
library(glmnet) 

library(tidyverse)
library(ggpubr) 
library(latex2exp)

library(foreach)
library(doParallel)


source(here("R", "utils.R"))
source(here("R", "methods.R"))
source(here("R", "plot.R"))


experiment <- "ROC_compare"

p <- 500
n <- 3*p
X_type <- "MCC_Block"
X_seed <- 2021

pi1 <- 50 / p
posit_type <- "random"
noise <- quote(rnorm(n))

target <- 0.3
target_at_alpha <- 0.2

sample_size <- 1

n_cores <- 14


X <- gene_X(X_type, n, p, X_seed)

mu1 <- BH_lm_calib(X, pi1, noise, posit_type, 1, side = "two", nreps = 200,
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
        return(c(FPP, FDP, TPP))
    })
    OLS_ROC <- data.frame(FPP = OLS_ROC[1, ], FDP = OLS_ROC[2, ],
                          TPP = OLS_ROC[3, ], stat = "OLS")
    
    LASSO_ROC <- sapply(sort(LASSO_beta), function(thr){
        selected <- which(LASSO_beta >= thr)
        FDP <- length(intersect(selected, which(H0))) / max(length(selected), 1)
        FPP <- length(intersect(selected, which(H0))) / ((1-pi1) * p)
        TPP <- length(intersect(selected, which(!H0))) / (pi1 * p)
        return(c(FPP, FDP, TPP))
    })
    LASSO_ROC <- data.frame(FPP = LASSO_ROC[1, ], FDP = LASSO_ROC[2, ],
                            TPP = LASSO_ROC[3, ], stat = "LASSO")
    
    return(list(OLS_ROC = OLS_ROC, LASSO_ROC = LASSO_ROC))
}

result <- results[[1]]
plot_data <- rbind(result$OLS_ROC, result$LASSO_ROC) %>%
    mutate(FPP = FPP + 1e-10 * TPP, FDP = FDP + 1e-10 * TPP)

method_colors <- c("dodgerblue3", "red")
method_names <- unname(TeX(c("$|\\hat{\\beta}_{OLS}|$", "$|\\hat{\\beta}_{LASSO}|$")))

plot <- ggplot(plot_data) +
    geom_line(aes(x = FPP, y = TPP, color = stat)) +
    scale_color_manual(values = method_colors, breaks = c("OLS", "LASSO"),
                       labels = method_names) +
    theme_bw() +
    theme(aspect.ratio = 1,
          # panel.grid = element_blank(),
          strip.text = element_text(size = 13),
          axis.title = element_text(size = 11),
          axis.text = element_text(size = 8),
          legend.position = "right",
          legend.title=element_text(size=9),
          legend.text=element_text(size=9)) +
    labs(color = "reject for large")

# figure <- ggarrange(plot1, plot2,
#                     # labels = c("A", "B", "C"),
#                     ncol = 1, nrow = 2,
#                     common.legend = T, legend = "right")

ggsave(filename = here("figs", paste0("ROC_compare.pdf")),
       plot, width = 4.5, height = 3.5)
# ggsave(filename = here("figs", paste0("ROC_compare.pdf")),
#        figure, width = 4, height = 6)


