
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

precompute_BonfBH_X <- function(X, knockoffs = knockoff::create.fixed){
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
    
    n <- NROW(X.pack$X)
    p <- NCOL(X.pack$X)
    df <- n-p
    
    y.org.nrow <- length(y)
    
    y <- matrix(y, nrow = 1)
    y_Pi_X <- X.pack$X_basis %*% t(y %*% X.pack$X_basis[1:y.org.nrow, ])
    y_Pi_X_res_norm2 <- sum(y^2) - sum(y_Pi_X^2)
    
    if(y.org.nrow < 2*p){
        if(randomize){
            y.extra <- rnorm(2*p-y.org.nrow,
                             sd = sqrt(y_Pi_X_res_norm2 / (X.pack$X.org.nrow - p)))
        } else{
            y.extra <- with_seed(0, rnorm(2*p-y.org.nrow,
                                          sd = sqrt(y_Pi_X_res_norm2 / (X.pack$X.org.nrow - p))))
        }
        y <- c(y, y.extra)
    }
    if(y.org.nrow < 2*p+1){
        y <- c(y, sqrt(y_Pi_X_res_norm2 / (X.pack$X.org.nrow - p)))
    }
    y_Pi_X_res_norm2 <- sum(y^2) - sum(y_Pi_X^2)
    
    y <- matrix(y, nrow = 1)
    
    
    vjy_obs <- c(y %*% X.pack$vj_mat)
    
    y_Pi_Xnoj <- sapply(1:p, function(j){
        y_Pi_X - X.pack$vj_mat[, j] * vjy_obs[j]
    })
    
    y_Pi_Xnoj_res_norm2 <- y_Pi_X_res_norm2 + vjy_obs^2
    sigmahat_XXk_res <- sqrt((y_Pi_X_res_norm2 - sum((y %*% X.pack$X_res_Xk_basis)^2)) / (n - 2*p))
    
    return(list(n = n, p = p, df = df, y = as.vector(y), vjy_obs = vjy_obs,
                y_Pi_X_res_norm2 = y_Pi_X_res_norm2, y_Pi_Xnoj_res_norm2 = y_Pi_Xnoj_res_norm2,
                sigmahat_XXk_res = sigmahat_XXk_res,
                y_Pi_Xnoj = y_Pi_Xnoj))
}


get_multi_method_list <- function(X, knockoffs, statistic){
    n <- NROW(X)
    p <- NCOL(X)
    
    Sigma <- solve(t(X) %*% X)
    X.pack <- process_X(X, knockoffs = knockoffs)
    if(n > 2*p){
        BonfBH_X <- precompute_BonfBH_X(X)
    } else{
        BonfBH_X <- NA
        warning('Input X has dimensions n < 2p+1. ',
                'Cannot use BonfBH.', immediate.=T)
    }
    
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
            y <- y.pack$y
            sigma_tilde <- y.pack$sigmahat_XXk_res
            
            kn_stats_obs <- statistic(X.pack$X, X.pack$X_kn, y, sigma_tilde = sigma_tilde)
            
            result <- kn.select(kn_stats_obs, alpha, selective = T, early_stop = F)
            sign_predict <- sign(c(matrix(y, nrow = 1) %*% (X.pack$X - X.pack$X_kn)))
            
            return(list(selected = result$selected, sign_predict = sign_predict))
        },
        BonBH = function(y, X, alpha){
            result <- BonfBH(X, y, alpha, BonfBH_X = BonfBH_X)
            return(list(selected = result$rejs))
        },
        cKnockoff = function(y, X, alpha){
            result <- cknockoff(X, y,
                                statistic = statistic,
                                alpha = alpha,
                                kappa = 1,
                                n_cores = 1,
                                X.pack = X.pack)
            return(list(selected = result$selected, sign_predict = result$sign_predict))
        },
        cKnockoff_L_0d8_R_ = function(y, X, alpha){
            result <- cknockoff(X, y,
                                statistic = statistic,
                                alpha = alpha,
                                kappa = 0.8,
                                n_cores = 1,
                                X.pack = X.pack)
            return(list(selected = result$selected, sign_predict = result$sign_predict))
        }
    )
    
    return(methods)
}

multi_method_color <- c("#984ea3", "dodgerblue3", "#333333", "#33a02c", "red", "orange1")
names(multi_method_color) <- c("BH", "dBH", "knockoff", "BonBH", "cKnockoff", "cKnockoff_L_0d8_R_")

multi_method_shape <- c(3, 4, 17, 23, 19, 15)
names(multi_method_shape) <- c("BH", "dBH", "knockoff", "BonBH", "cKnockoff", "cKnockoff_L_0d8_R_")


get_kn_method_list <- function(X, knockoffs, statistic){
    X.pack <- process_X(X, knockoffs = knockoffs)
    knockoffs_gene <- function(X){return(X.pack$X_kn)}
    
    methods <- list(
        kn_D_lambdasmax = function(y, X, alpha){
            y.pack <- process_y(X.pack, y)
            y <- y.pack$y
            
            result <- knockoff.filter(X.pack$X, y, knockoffs = knockoffs_gene,
                                      statistic = stat.glmnet_lambdasmax, fdr = alpha)
            sign_predict <- sign(c(matrix(y, nrow = 1) %*% (X.pack$X - X.pack$X_kn)))
            
            return(list(selected = result$selected, sign_predict = sign_predict))
        },
        kn_D_lambdasmax_lm = function(y, X, alpha){
            y.pack <- process_y(X.pack, y)
            y <- y.pack$y
            sigma_tilde <- y.pack$sigmahat_XXk_res
            
            kn_stats_obs <- stat.glmnet_lambdasmax_lm(X.pack$X, X.pack$X_kn, y, sigma_tilde = sigma_tilde)
            
            result <- kn.select(kn_stats_obs, alpha, selective = T, early_stop = F)
            sign_predict <- sign(c(matrix(y, nrow = 1) %*% (X.pack$X - X.pack$X_kn)))
            
            return(list(selected = result$selected, sign_predict = sign_predict))
        },
        kn_D_coefdiff_lm = function(y, X, alpha){
            y.pack <- process_y(X.pack, y)
            y <- y.pack$y
            sigma_tilde <- y.pack$sigmahat_XXk_res
            
            kn_stats_obs <- stat.glmnet_coefdiff_lm(X.pack$X, X.pack$X_kn, y, sigma_tilde = sigma_tilde)
            
            result <- kn.select(kn_stats_obs, alpha, selective = T, early_stop = F)
            sign_predict <- sign(c(matrix(y, nrow = 1) %*% (X.pack$X - X.pack$X_kn)))
            
            return(list(selected = result$selected, sign_predict = sign_predict))
        }
    )
    
    return(methods)
}


kn_method_color <- c("#1b9e77", "#d95f02", "#7570b3")
names(kn_method_color) <- c("kn_D_lambdasmax", "kn_D_lambdasmax_lm", "kn_D_coefdiff_lm")

kn_method_shape <- c(4, 1, 2)
names(kn_method_shape) <- c("kn_D_lambdasmax", "kn_D_lambdasmax_lm", "kn_D_coefdiff_lm")

