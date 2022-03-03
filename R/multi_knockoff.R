# borrowed from Lihua https://github.com/lihualei71/mkn/blob/master/R/solveS.R
mineig <- function(A){
    min(eigen(A, symmetric = TRUE)$values)
}

# borrowed from Lihua https://github.com/lihualei71/mkn/blob/master/R/solveS.R
mkn_solve_sdp <- function(Sigma, k, gaptol = 1e-06,
                          maxit = 1000, psdtol = 1e-09){
    stopifnot(isSymmetric(Sigma))
    ratio <- (k + 1) / k
    G <- stats::cov2cor(Sigma)
    p <- dim(G)[1]
    if (mineig(G) < psdtol) {
        stop("The covariance matrix is not positive-definite: cannot solve SDP", 
             immediate. = T)
    }
    Cl1 <- rep(0, p)
    Al1 <- -Matrix::Diagonal(p)
    Cl2 <- rep(1, p)
    Al2 <- Matrix::Diagonal(p)
    d_As <- c(diag(p))
    As <- Matrix::Diagonal(length(d_As), x = d_As)
    As <- As[which(Matrix::rowSums(As) > 0), ]
    Cs <- c(ratio * G)
    A <- cbind(Al1, Al2, As)
    C <- matrix(c(Cl1, Cl2, Cs), 1)
    K <- NULL
    K$s <- p
    K$l <- 2 * p
    b <- rep(1, p)
    OPTIONS <- NULL
    OPTIONS$gaptol <- gaptol
    OPTIONS$maxit <- maxit
    OPTIONS$logsummary <- 0
    OPTIONS$outputstats <- 0
    OPTIONS$print <- 0
    sol <- Rdsdp::dsdp(A, b, C, K, OPTIONS)
    if (!identical(sol$STATS$stype, "PDFeasible")) {
        warning("The SDP solver returned a non-feasible solution")
    }
    s <- sol$y
    s[s < 0] <- 0
    s[s > 1] <- 1
    psd <- 0
    s_eps <- 1e-08
    while (mineig(ratio * G - diag(s * (1 - s_eps), length(s))) < psdtol) {
        s_eps <- s_eps * 10
    }
    s <- s * (1 - s_eps)
    if (max(s) == 0) {
        warning("In creation of SDP knockoffs, procedure failed. Knockoffs will have no power.", 
                immediate. = T)
    }
    return(s * diag(Sigma))
}

create.mkn <- function(X, k){
    n <- NROW(X)
    p <- NCOL(X)
    
    X <- scale(X, center = FALSE, scale = sqrt(colSums(X^2)))
    
    Sigma <- t(X) %*% X
    S <- mkn_solve_sdp(Sigma, k)
    
    diag_mat <- diag(nrow = k+1)
    offdiag_mat <- matrix(1, nrow = k+1, ncol = k+1) - diag_mat
    Gram <- diag_mat %x% Sigma + offdiag_mat %x% (Sigma - diag(S))
    
    svd_res <- svd(Gram)
    X0 <- svd_res$u %*% diag(sqrt(svd_res$d)) %*% t(svd_res$u)
    
    R0 <- qr.R(qr(X0))
    QR <- qr(X)
    Q <- qr.Q(QR, complete = T)
    R0 <- diag(sign(diag(qr.R(QR))) * sign(diag(R0))) %*% R0
    
    X1 <- Q[, 1:((k+1)*p)] %*% R0
    
    return(list(X = X, Xk = X1[, (p+1):((k+1)*p)],
                res_basis = Q[, ((k+1)*p+1):n]))
}


score.glmnet_coefdiff_lm <- function(X, X_k, y, sigma_tilde) {
    n <- NROW(X)
    p <- NCOL(X)
    k <- NCOL(X_k) / p
    orig <- 1:p
    
    lambda <- 2 * sigma_tilde / n
    
    # Randomly swap columns of X and Xk
    swap_index_mat <- sapply(1:p, function(feature_id){
        sample(0:k) * p + feature_id
    })
    swap_index <- c(t(swap_index_mat))
    X_full <- cbind(X, X_k)[, swap_index]
    
    # Compute statistics
    fit <- glmnet::glmnet(X_full, y, lambda = lambda,
                          intercept = T, standardize = F,
                          standardize.response = F, family = "gaussian")
    
    Z <- abs(c(as.matrix(fit$beta)))
    
    y_res <- as.matrix(y - (X_full %*% fit$beta + rep(fit$a0, each = n)))
    LM <- c(abs(t(y_res) %*% X_full))
    LM_add_on <- rep(lambda * n, length(fit$beta))
    not_in_model <- which(abs(fit$beta) < 1e-10)
    LM_add_on[not_in_model] <- LM[not_in_model]
    
    scores <- abs(fit$beta) + LM_add_on
    
    order_stats <- sapply(1:p, function(feature_id){
        max(scores[(0:k) * p + feature_id])
    })
    org_ranks <- sapply(1:p, function(feature_id){
        score_aug <- scores[(0:k) * p + feature_id]
        org_ind <- as.integer((c(swap_index_mat[, feature_id]) - feature_id) / p)
        org_ind <- which(org_ind == 0)
        rank_org <- rank(-score_aug)[org_ind]
        return(rank_org)
    })
    
    return(list(order_stats = order_stats, org_ranks = org_ranks))
}

mkn.select <- function(order_stats, org_ranks, k, alpha){
    p <- length(order_stats)
    s <- 1
    lambda <- ceiling((k+1)/2)
    
    order_desc_ind <- order(order_stats, decreasing = T)
    rank_ordered <- org_ranks[order_desc_ind]
    
    # fdp <- sapply(abs(W_desc), function(t){
    #     (1 + sum(kn_statistics <= -t)) / max(1, sum(kn_statistics >= t))
    # })
    fdp <- rep(NA, p)
    for(iter in 1:p){
        pos_stat_num <- sum(rank_ordered[1:iter] <= s)
        neg_stat_num <- sum(rank_ordered[1:iter] > lambda)
        fdp[iter] <- (1 + neg_stat_num) / max(1, pos_stat_num) * s / (k+1 - lambda)
    }
    ok <- which(fdp <= alpha)
    k_hat <- ifelse(length(ok) > 0, max(ok), 0)
    
    if(k_hat > 0){
        selected <- sort(order_desc_ind[which(rank_ordered[1:k_hat] <= s)])
    } else{
        selected <- NULL
    }
    fdp_est <- ifelse(k_hat > 0, fdp[k_hat], 1)
    
    results <- list(selected = selected, fdp_est = fdp_est)
}



mkn.filter <- function(X, y, k, alpha){
    
}