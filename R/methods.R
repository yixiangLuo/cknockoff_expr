source(here("R", "multi_knockoff.R"))


# weighted BH method
BH_weighted <- function(pvals, alpha,
                        weights = rep(1, length(pvals)) / length(pvals)){
    n <- length(pvals)

    adjust_pvals <- sort(pvals / weights) / (1:n)
    nrejs <- max(0, which(adjust_pvals <= alpha))

    rejs <- which(pvals <= nrejs * alpha * weights)

    return(list(nrejs = nrejs, rejs = rejs))
}

# compute the p-values of t-statistics
pvals_t <- function(tvals, df, side = "two"){
    if (side == "right"){
        pvals <- pt(tvals, df = df, lower.tail = FALSE)
    } else if (side == "left"){
        pvals <- pt(tvals, df = df, lower.tail = TRUE)
    } else if (side == "two"){
        pvals <- 2 * pt(abs(tvals), df = df, lower.tail = FALSE)
    }

    return(pvals)
}

# compute the t-vals of a linear regression test problem
lm_to_t <- function(y, X, Sigma = NULL){
    n <- NROW(X)
    p <- NCOL(X)

    if(is.null(Sigma)){
        Sigma <- solve(t(X) %*% X)
    }
    Xy <- t(X) %*% y
    df <- n - p

    zvals <- Sigma %*% Xy
    sigmahat <- as.numeric(sqrt((sum(y^2) - t(Xy) %*% zvals) / df))
    tvals <- zvals / sqrt(diag(Sigma)) / sigmahat

    return(list(tvals = tvals, df = df))
}

# apply BH method to a linear regression test problem
BH_lm <- function(y, X, side = "two", alpha,
                  weights = rep(1, NCOL(X)) / NCOL(X),
                  Sigma = NULL){

    t_result <- lm_to_t(y, X, Sigma)
    pvals <- pvals_t(t_result$tvals, t_result$df, side = "two")

    BH_result <-  BH_weighted(pvals, alpha, weights)

    return(BH_result)
}

precompute_BonfBH_X <- function(X, knockoffs = ckn.create.fixed){
    knockoff <- knockoffs(X)
    X <- knockoff$X
    X_kn <- knockoff$Xk
    # X_full <- cbind(X, X_kn)
    # X_full_H <- X_full %*% solve(t(X_full) %*% X_full) %*% t(X_full)

    kn_diff_D <- pmax(diag(t(X) %*% X - t(X_kn) %*% X), 0)

    A_inv <- solve(2 * t(X) %*% X - diag(kn_diff_D))

    return(list(X = X, X_kn = X_kn, kn_diff_D = kn_diff_D, A_inv = A_inv))
}

BonfBH <- function(X, y, alpha, BonfBH_X = NULL){
    if(is.null(BonfBH_X)){
        BonfBH_X <- precompute_BonfBH_X(X)
    }
    X <- BonfBH_X$X

    n <- NROW(X)
    p <- NCOL(X)
    df <- n - 2*p

    hbeta_1 <- BonfBH_X$A_inv %*% t(X + BonfBH_X$X_kn) %*% y
    hbeta_2 <- t(X - BonfBH_X$X_kn) %*% y / BonfBH_X$kn_diff_D

    y_proj_X_full <- lm(y ~ cbind(X, BonfBH_X$X_kn) + 0)$fitted
    sigmahat <- as.numeric(sqrt((sum(y^2) - sum(y_proj_X_full^2)) / df))

    tvals1 <- hbeta_1 / sqrt(2 * diag(BonfBH_X$A_inv)) / sigmahat
    tvals2 <- hbeta_2 / sqrt(2 / BonfBH_X$kn_diff_D) / sigmahat

    pvals1 <- pvals_t(tvals1, df, side = "two")
    pvals1[is.na(pvals1)] <- 1
    pvals2 <- pvals_t(tvals2, df, side = "two")
    pvals2[is.na(pvals2)] <- 1

    pi0_hat <- 2 * (sum(pvals2 > 0.5) + 1) / p
    p_star <- pvals2 * pi0_hat
    p_star[pvals1 > sqrt(alpha)] <- 1

    result <- BH_weighted(p_star, sqrt(alpha))

    return(result)
}




process_y <- function(X.pack, y, randomize = F){

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
    
    sigmahat_XXk_res <- sqrt((RSS_XXk) / df_XXk)

    return(list(y = as.vector(y), y.org = as.vector(y.org),
                sigmahat_XXk_res = sigmahat_XXk_res))
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
    # if(early_stop == 2){
    #   ok <- which((fdp <= 1) & (cumsum(W_desc > 0) < 1/alpha))
    #   k_hat <- max(k_hat, ifelse(length(ok) > 0, max(ok), 0))
    # }
    
    if(k_hat > 0){
        selected <- sort(W_desc_ind[which(W_desc[1:k_hat] > 0)])
    } else{
        selected <- NULL
    }
    fdp_est <- ifelse(k_hat > 0, fdp[k_hat], 1)
    W_k_hat <- abs(W_desc[max(k_hat, 1)])
    
    results <- list(selected = selected, fdp_est = fdp_est,  W_k_hat = W_k_hat)
}


get_multi_method_list <- function(X, knockoffs, statistic){
    n <- NROW(X)
    p <- NCOL(X)

    Sigma <- solve(t(X) %*% X)
    X.pack <- process_X(X, knockoffs = knockoffs, intercept = F)
    if(n > 2*p){
        BonfBH_X <- precompute_BonfBH_X(X)
    } else{
        BonfBH_X <- NA
        warning('Input X has dimensions n < 2p+1. ',
                'Cannot use BonfBH.', immediate.=T)
    }
    # mkn_k <- 5
    # mkn_X.pack <- create.mkn(X, mkn_k)

    methods <- list(
        BH = function(y, X, alpha){
            # result <- BH_lm(y, X, side = "two", alpha, Sigma = Sigma)

            t_result <- lm_to_t(y, X, Sigma)
            pvals <- pvals_t(t_result$tvals, t_result$df, side = "two")
            result <-  BH_weighted(pvals, alpha)

            sign_predict <- rep(0, NCOL(X))
            sign_predict[result$rejs] <- sign(t_result$tvals[result$rejs])

            return(list(selected = result$rejs, sign_predict = sign_predict))
        },
        dBH = function(y, X, alpha){
            result <- dBH_lm(y, X, intercept = FALSE, side = "two",
                             alpha = alpha, gamma = 0.9, niter = 1,
                             avals_type = "BH", qcap = 2)

            t_result <- lm_to_t(y, X, Sigma)
            sign_predict <- rep(0, NCOL(X))
            sign_predict[result$rejs] <- sign(t_result$tvals[result$rejs])

            return(list(selected = result$rejs, sign_predict = sign_predict))
        },
        knockoff = function(y, X, alpha){
            y.pack <- process_y(X.pack, y)
            sigma_tilde <- y.pack$sigmahat_XXk_res
            
            kn_stats_obs <- statistic(X.pack$X.org, X.pack$X_kn.org,
                                      y.pack$y.org, sigma_tilde = sigma_tilde)

            result <- kn.select(kn_stats_obs, alpha, selective = T, early_stop = F)
            sign_predict <- sign(c(matrix(y.pack$y.org, nrow = 1) %*% 
                                       (X.pack$X.org - X.pack$X_kn.org)))

            return(list(selected = result$selected, sign_predict = sign_predict))
        },
        # mKnockoff = function(y, X, alpha){
        #     sigma_tilde <- sqrt((sum((matrix(y, nrow = 1) %*% mkn_X.pack$res_basis)^2)) / (n - (mkn_k+1)*p))
        #     
        #     scores <- score.glmnet_coefdiff_lm(mkn_X.pack$X, mkn_X.pack$Xk, y, sigma_tilde)
        #     
        #     result <- mkn.select(scores$order_stats, scores$org_ranks, mkn_k, alpha)
        #     sign_predict <- rep(NA, NCOL(X))
        #     
        #     return(list(selected = result$selected, sign_predict = sign_predict))
        # },
        BonBH = function(y, X, alpha){
            result <- BonfBH(X, y, alpha, BonfBH_X = BonfBH_X)
            return(list(selected = result$rejs))
        },
        cKnockoff = function(y, X, alpha){
            result <- cknockoff(X, y,
                                intercept = F,
                                statistic = statistic,
                                alpha = alpha,
                                n_cores = 1,
                                X.pack = X.pack,
                                Rstar_refine = F)
            return(list(selected = result$selected, sign_predict = result$sign_predict))
        },
        cKnockoff_STAR = function(y, X, alpha){
            result <- cknockoff(X, y,
                                intercept = F,
                                statistic = statistic,
                                alpha = alpha,
                                n_cores = 1,
                                X.pack = X.pack,
                                Rstar_refine = T)
            return(list(selected = result$selected, sign_predict = result$sign_predict))
        }
    )

    return(methods)
}

multi_method_color <- c("#984ea3", "dodgerblue3", "#333333", "#006d2c", "#33a02c", "red", "orange1")
names(multi_method_color) <- c("BH", "dBH", "knockoff", "mKnockoff", "BonBH", "cKnockoff", "cKnockoff_STAR")

multi_method_shape <- c(3, 4, 17, 6, 23, 19, 15)
names(multi_method_shape) <- c("BH", "dBH", "knockoff", "mKnockoff", "BonBH", "cKnockoff", "cKnockoff_STAR")


get_kn_method_list <- function(X, knockoffs, statistic){
    X.pack <- process_X(X, knockoffs = knockoffs, intercept = F)
    knockoffs_gene <- function(X){return(X.pack$X_kn)}

    methods <- list(
        kn_D_lambdasmax = function(y, X, alpha){
            y.pack <- process_y(X.pack, y)

            result <- knockoff.filter(X.pack$X.org, y.pack$y.org, knockoffs = knockoffs_gene,
                                      statistic = stat.glmnet_lambdasmax, fdr = alpha)
            sign_predict <- sign(c(matrix(y.pack$y.org, nrow = 1) %*% 
                                       (X.pack$X.org - X.pack$X_kn.org)))

            return(list(selected = result$selected, sign_predict = sign_predict))
        },
        kn_D_lambdasmax_lm = function(y, X, alpha){
            y.pack <- process_y(X.pack, y)
            sigma_tilde <- y.pack$sigmahat_XXk_res

            kn_stats_obs <- stat.glmnet_lambdasmax_lm(X.pack$X.org, X.pack$X_kn.org,
                                                      y.pack$y.org, sigma_tilde = sigma_tilde)

            result <- kn.select(kn_stats_obs, alpha, selective = T, early_stop = F)
            sign_predict <- sign(c(matrix(y.pack$y.org, nrow = 1) %*% 
                                       (X.pack$X.org - X.pack$X_kn.org)))

            return(list(selected = result$selected, sign_predict = sign_predict))
        },
        kn_D_coefdiff_lm = function(y, X, alpha){
            y.pack <- process_y(X.pack, y)
            sigma_tilde <- y.pack$sigmahat_XXk_res

            kn_stats_obs <- stat.glmnet_coefdiff_lm(X.pack$X.org, X.pack$X_kn.org,
                                                    y.pack$y.org, sigma_tilde = sigma_tilde)

            result <- kn.select(kn_stats_obs, alpha, selective = T, early_stop = F)
            sign_predict <- sign(c(matrix(y.pack$y.org, nrow = 1) %*% 
                                       (X.pack$X.org - X.pack$X_kn.org)))

            return(list(selected = result$selected, sign_predict = sign_predict))
        }
    )

    return(methods)
}


kn_method_color <- c("#1b9e77", "#d95f02", "#7570b3")
names(kn_method_color) <- c("kn_D_lambdasmax", "kn_D_lambdasmax_lm", "kn_D_coefdiff_lm")

kn_method_shape <- c(4, 1, 2)
names(kn_method_shape) <- c("kn_D_lambdasmax", "kn_D_lambdasmax_lm", "kn_D_coefdiff_lm")

