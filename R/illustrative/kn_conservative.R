library(here)
library(abind) 

library(tidyverse)

library(knockoff)
library(cknockoff)

library(foreach)
library(doParallel)


source(here("R", "utils.R"))
source(here("R", "methods.R"))
source(here("R", "plot.R"))


experiment <- "kn_conservative"

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

knockoffs <- create.fixed
statistic <- stat.glmnet_coefdiff_lm

# alphas <- seq(from = 0.005, to = 0.3, length.out = 100)
alphas <- seq(from = 0.05, to = 0.2, by = 0.05)



X <- gene_X(X_type, n, p, X_seed)
X.pack <- process_X(X, knockoffs = knockoffs, intercept = F)

mu1 <- BH_lm_calib(X, pi1, noise, posit_type, 1, side = "two", nreps = 200,
                   alpha = target_at_alpha, target = target, n_cores = n_cores)
beta <- genmu(p, pi1, mu1, posit_type, 1)

H0 <- beta == 0


registerDoParallel(n_cores)



transform_y <- function(X.pack, y, randomize = F){
    n <- length(y)
    p <- NCOL(X.pack$X)
    
    # compute RSS and degree of freedom in y ~ X
    RSS_X <- sum(y^2) - sum((matrix(y, nrow=1) %*% X.pack$XXk.org.basis[1:n, 1:p])^2)
    df_X <- n - p
    
    # augment y if needed
    if(n < 2*p){
        if(randomize){
            y.extra <- rnorm(2*p-n, sd = sqrt(RSS_X / (n - p)))
        } else{
            y.extra <- with_seed(0, rnorm(2*p-n, sd = sqrt(RSS_X / (n - p))))
        }
        y <- c(y, y.extra)
    }
    
    # record the original y before rotation
    y.org <- y
    
    # rotate y in the same way as we did for [X Xk]
    y.norm2 <- sum(y^2)
    y <- as.vector(matrix(y, nrow=1) %*% X.pack$XXk.org.basis)
    
    # the component of y not in the 2p-dim subspace (residue of y~[X Xk])
    # is recorded separately as RSS_XXk and df_XXk
    if(n < 2*p+1){
        RSS_XXk <- RSS_X
        df_XXk <- df_X
    } else{
        RSS_XXk <- y.norm2 - sum(y^2)
        df_XXk <- n - 2*p
    }
    
    y.data <- list(y = y, y.org = y.org,
                   RSS_X = RSS_X, RSS_XXk = RSS_XXk,
                   df_X = df_X, df_XXk = df_XXk)
    return(y.data)
    
}

kn.select <- function(kn_statistics, alpha,
                      selective = T, early_stop = F){
    p <- length(kn_statistics)
    W_desc_ind <- order(abs(kn_statistics), decreasing = T)
    W_desc <- kn_statistics[W_desc_ind]
    
    # fdp <- sapply(abs(W_desc), function(t){
    #     (1 + sum(kn_statistics <= -t)) / max(1, sum(kn_statistics >= t))
    # })
    fdp <- rep(NA, p)
    neg_stat_num <- 0
    pos_stat_num <- 0
    for(k in 1:p){
        pos_stat_num <- pos_stat_num + (W_desc[k] > 0)
        neg_stat_num <- neg_stat_num + (W_desc[k] < 0)
        fdp[k] <- (1 + neg_stat_num) / max(1, pos_stat_num)
    }
    ok <- which(fdp <= alpha)
    k_hat <- ifelse(length(ok) > 0, max(ok), 0)
    
    if(early_stop){
        ok <- which(cumsum(W_desc > 0) < 1/alpha)
        k_hat <- max(k_hat, ifelse(length(ok) > 0, max(ok), 0))
    }
    
    if(k_hat > 0){
        selected <- sort(W_desc_ind[which(W_desc[1:k_hat] > 0)])
        candidates <- sort(W_desc_ind[1:k_hat])
    } else{
        selected <- NULL
        candidates <- NULL
    }
    fdp_est <- ifelse(k_hat > 0, fdp[k_hat], 1)
    W_k_hat <- abs(W_desc[max(k_hat, 1)])
    
    results <- list(selected = selected, candidates = candidates,
                    fdp_est = fdp_est,  W_k_hat = W_k_hat)
}




print(Sys.time())

results <- foreach(iter = 1:sample_size) %dopar% {
    
    set.seed(iter)
    
    y <- X %*% beta + eval(noise)
    
    y.data <- transform_y(X.pack, y)
    
    kn.stats <- statistic(X.pack$X, X.pack$X_kn, y.data$y,
                          sigma_tilde = sqrt(y.data$RSS_XXk / y.data$df_XXk))
    
    result_y <- sapply(alphas, function(alpha){
        kn.result <- kn.select(kn.stats, alpha, selective = T, early_stop = F)
        kn.early.result <- kn.select(kn.stats, alpha, selective = T, early_stop = T)
        
        FDP <- length(intersect(kn.result$selected, which(H0))) / max(length(kn.result$selected), 1)
        TDP <- length(intersect(kn.result$selected, which(!H0))) / max(length(kn.result$selected), 1)
        
        b0.F <- FDP * alpha / kn.result$fdp_est
        b0.T <- TDP * alpha / kn.result$fdp_est
        
        FDP_early <- length(intersect(kn.early.result$selected, which(H0))) / max(length(kn.early.result$selected), 1)
        TDP_early <- length(intersect(kn.early.result$selected, which(!H0))) / max(length(kn.early.result$selected), 1)
        
        b.F <- FDP_early * alpha / kn.early.result$fdp_est
        b.T <- TDP_early * alpha / kn.early.result$fdp_est
        
        M_tau0 <- sum((kn.stats[kn.result$candidates] > 0) & H0[kn.result$candidates]) / 
            (1 + sum((kn.stats[kn.result$candidates] < 0) & H0[kn.result$candidates]))
        M_tau <- sum((kn.stats[kn.early.result$candidates] > 0) & H0[kn.early.result$candidates]) / 
            (1 + sum((kn.stats[kn.early.result$candidates] < 0) & H0[kn.early.result$candidates]))
        
        M0 <- sum((kn.stats > 0) & H0) / (1 + sum((kn.stats < 0) & H0))
        
        result <- c(M0, M_tau, b.F / alpha, M_tau0, b0.F / alpha, FDP / alpha,
                    b.T / alpha, b0.T / alpha, TDP / alpha)
        names(result) <- c("M_0", "M_tau", "b_F", "M_tau0", "b0_F", "FDR",
                           "b_T", "b0_T", "TDR")
        
        return(result)
    })
    
    return(result_y)
}

print(Sys.time())

results <- abind(results, along=3)
results <- apply(results, c(1, 2), mean)
results <- rbind(results, alpha = alphas)

save(results, file = here("data", paste0(experiment, ".Rdata")))


draw_kn_conservative(experiment)





