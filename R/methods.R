library(here)
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
    # pvals1[is.na(pvals1)] <- 1
    pvals2 <- pvals_t(tvals2, df, side = "two")
    # pvals2[is.na(pvals2)] <- 1

    pi0_hat <- 2 * (sum(pvals2 > 0.5) + 1) / p
    p_star <- pvals2 * pi0_hat
    p_star[pvals1 > sqrt(alpha)] <- 1

    result <- BH_weighted(p_star, sqrt(alpha))

    return(result)
}


get_multi_method_list <- function(X, knockoffs, statistic, method_names, Xcov.true = NA){
    
    methods <- list()
    
    # preprocess if needed
    if("BH" %in% method_names | "dBH" %in% method_names){
        Sigma <- solve(t(X) %*% X)
    }
    if("knockoff" %in% method_names | "cKnockoff" %in% method_names | 
       "cKnockoff_STAR" %in% method_names){
        X.pack <- cknockoff::process_X(X, knockoffs = knockoffs, intercept = F)
    }
    if("BonBH" %in% method_names){
        BonfBH_X <- precompute_BonfBH_X(X)
    }
    if("mKnockoff" %in% method_names){
        mkn_k <- 5
        n <- NROW(X)
        p <- NCOL(X)
        
        if(n > mkn_k * p){
            mkn_X.pack <- create.mkn(X, mkn_k)
        }
    }
    
    # create method functions
    if("BH" %in% method_names){
        methods$BH <- function(y, X, alpha){
            # result <- BH_lm(y, X, side = "two", alpha = alpha, Sigma = Sigma)
            
            t_result <- lm_to_t(y, X, Sigma)
            pvals <- pvals_t(t_result$tvals, t_result$df, side = "two")
            result <-  BH_weighted(pvals, alpha)
            
            sign_predict <- rep(0, NCOL(X))
            sign_predict[result$rejs] <- sign(t_result$tvals[result$rejs])
            
            return(list(selected = result$rejs, sign_predict = sign_predict))
        }
    }
    if("dBH" %in% method_names){
        methods$dBH <- function(y, X, alpha){
            result <- dBH_lm(y, X, intercept = FALSE, side = "two",
                             alpha = alpha, gamma = 0.9, niter = 1,
                             avals_type = "BH", qcap = 2)
            
            t_result <- lm_to_t(y, X, Sigma)
            sign_predict <- rep(0, NCOL(X))
            sign_predict[result$rejs] <- sign(t_result$tvals[result$rejs])
            
            return(list(selected = result$rejs, sign_predict = sign_predict))
        }
    }
    if("knockoff" %in% method_names){
        methods$knockoff <- function(y, X, alpha){
            y.data <- cknockoff:::transform_y(X.pack, y, intercept = F)
            sigma_tilde <- sqrt(y.data$RSS_XXk / y.data$df_XXk)
            
            kn_stats_obs <- statistic(X.pack$X.org, X.pack$X_kn.org,
                                      y.data$y.org, sigma_tilde = sigma_tilde)
            
            result <- cknockoff:::kn.select(kn_stats_obs, alpha, selective = T, early_stop = F)
            sign_predict <- sign(c(matrix(y.data$y.org, nrow = 1) %*% 
                                       (X.pack$X.org - X.pack$X_kn.org)))
            
            return(list(selected = result$selected, sign_predict = sign_predict))
        }
    }
    if("mKnockoff" %in% method_names){
        if(n > mkn_k * p){
            methods$mKnockoff <- function(y, X, alpha){
                sigma_tilde <- sqrt((sum((matrix(y, nrow = 1) %*% mkn_X.pack$res_basis)^2)) / (n - (mkn_k+1)*p))
                
                scores <- score.glmnet_coefdiff_lm(mkn_X.pack$X, mkn_X.pack$Xk, y, sigma_tilde)
                
                result <- mkn.select(scores$order_stats, scores$org_ranks, mkn_k, alpha)
                sign_predict <- rep(NA, NCOL(X))
                
                return(list(selected = result$selected, sign_predict = sign_predict))
            }
        } else{
            methods$mKnockoff <- function(y, X, alpha){
                return(list(selected = NULL, sign_predict = NULL))
            }
        }
        
    }
    if("BonBH" %in% method_names){
        methods$BonBH <- function(y, X, alpha){
            result <- BonfBH(X, y, alpha, BonfBH_X = BonfBH_X)
            return(list(selected = result$rejs))
        }
    }
    if("cKnockoff" %in% method_names){
        methods$cKnockoff <- function(y, X, alpha){
            result <- cknockoff(X, y,
                                intercept = F,
                                statistic = statistic,
                                alpha = alpha,
                                n_cores = 1,
                                X.pack = X.pack,
                                Rstar_refine = F)
            return(list(selected = result$selected, sign_predict = result$sign_predict))
        }
    }
    if("cKnockoff_STAR" %in% method_names){
        methods$cKnockoff_STAR <- function(y, X, alpha){
            result <- cknockoff(X, y,
                                intercept = F,
                                statistic = statistic,
                                alpha = alpha,
                                n_cores = 1,
                                X.pack = X.pack,
                                Rstar_refine = T)
            return(list(selected = result$selected, sign_predict = result$sign_predict))
        }
    }
    if("cLasso" %in% method_names){
        methods$cLasso <- function(y, X, alpha){
            result <- cLasso(X, y, alpha)
            return(list(selected = result$selected))
        }
    }
    
    if("knockoff.MX" %in% method_names){
        methods$knockoff.MX <- function(y, X, alpha){
            Sigma.inv <- solve(Xcov.true)
            s <- knockoff:::create.solve_sdp(Xcov.true)
            s[s <= 1e-5] <- 0
            
            Xk_mumat <- diag(p) - Sigma.inv %*% diag(s)
            mu_kn <- X %*% Xk_mumat
            cov_kn <- 2 * diag(s) - diag(s) %*% Sigma.inv %*% diag(s)
            cov_kn_half <- chol(cov_kn)
            
            X_kn <<- ckn.modelX:::rnorm_mult(mu_kn, cov_kn_half)
            
            if("sigma_tilde" %in% names(formals(statistic))){
                kn_stats_obs <- statistic(X, X_kn, y, sigma_tilde = 1)
            } else{
                kn_stats_obs <- statistic(X, X_kn, y)
            }
            
            result <- ckn.modelX:::kn.select(kn_stats_obs, alpha, selective = T, early_stop = F)
            sign_predict <- rep(NA, p)
            
            return(list(selected = result$selected, sign_predict = sign_predict))
        }
    }
    if("cKnockoff.MX" %in% method_names){
        methods$cKnockoff.MX <- function(y, X, alpha){
            result <- ckn.modelX::cknockoff.MX(X, y, Xcov.true, X_kn,
                                               statistic = statistic,
                                               alpha = alpha,
                                               family = "gaussian",   # gaussian, binomial
                                               n_cores = 1)
            return(list(selected = result$selected, sign_predict = rep(NA, p)))
        }
    }

    return(methods)
}

multi_method_color <- c("#984ea3", "dodgerblue3", "#333333", "#006d2c", "#33a02c", "red", "orange1", "orange1", "#333333", "red")
names(multi_method_color) <- c("BH", "dBH", "knockoff", "mKnockoff", "BonBH", "cKnockoff", "cKnockoff_STAR", "cLasso", "knockoff.MX", "cKnockoff.MX")

multi_method_shape <- c(3, 4, 17, 6, 23, 19, 15, 19, 17, 19)
names(multi_method_shape) <- c("BH", "dBH", "knockoff", "mKnockoff", "BonBH", "cKnockoff", "cKnockoff_STAR", "cLasso", "knockoff.MX", "cKnockoff.MX")


get_kn_method_list <- function(X, knockoffs, statistic, method_names, Xcov.true = NA){
    X.pack <- process_X(X, knockoffs = knockoffs, intercept = F)
    knockoffs_gene <- function(X){return(X.pack$X_kn.org)}
    
    lambda_scale <- 1

    methods <- list(
        LSM = function(y, X, alpha){
            y.data <- cknockoff:::transform_y(X.pack, y, intercept = F)

            result <- knockoff.filter(X.pack$X.org, y.data$y.org, knockoffs = knockoffs_gene,
                                      statistic = stat.glmnet_lambdasmax, fdr = alpha)
            sign_predict <- sign(c(matrix(y.data$y.org, nrow = 1) %*% 
                                       (X.pack$X.org - X.pack$X_kn.org)))

            return(list(selected = result$selected, sign_predict = sign_predict))
        },
        C_D_LSM = function(y, X, alpha){
            y.data <- cknockoff:::transform_y(X.pack, y, intercept = F)
            sigma_tilde <- sqrt(y.data$RSS_XXk / y.data$df_XXk) * lambda_scale

            kn_stats_obs <- stat.glmnet_lambdasmax_coarse(X.pack$X.org, X.pack$X_kn.org,
                                                          c(y.data$y.org), sigma_tilde = sigma_tilde)

            result <- cknockoff:::kn.select(kn_stats_obs, alpha, selective = T, early_stop = F)
            sign_predict <- sign(c(matrix(y.data$y.org, nrow = 1) %*% 
                                       (X.pack$X.org - X.pack$X_kn.org)))

            return(list(selected = result$selected, sign_predict = sign_predict))
        },
        LCD_D_T = function(y, X, alpha){
            y.data <- cknockoff:::transform_y(X.pack, y, intercept = F)
            sigma_tilde <- sqrt(y.data$RSS_XXk / y.data$df_XXk) * lambda_scale

            kn_stats_obs <- stat.glmnet_coefdiff_tiebreak(X.pack$X.org, X.pack$X_kn.org,
                                                          y.data$y.org, sigma_tilde = sigma_tilde)

            result <- cknockoff:::kn.select(kn_stats_obs, alpha, selective = T, early_stop = F)
            sign_predict <- sign(c(matrix(y.data$y.org, nrow = 1) %*% 
                                       (X.pack$X.org - X.pack$X_kn.org)))

            return(list(selected = result$selected, sign_predict = sign_predict))
        }
    )

    return(methods[method_names])
}


kn_method_color <- c("#1b9e77", "#d95f02", "#7570b3")
names(kn_method_color) <- c("LSM", "C_D_LSM", "LCD_D_T")

kn_method_shape <- c(4, 1, 2)
names(kn_method_shape) <- c("LSM", "C_D_LSM", "LCD_D_T")

