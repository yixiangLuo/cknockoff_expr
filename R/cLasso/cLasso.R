library(here)
source(here("R", "cLasso", "homotopy_lasso.R"))
source(here("R", "cLasso", "lasso_hFDR.R"))
source(here("R", "cLasso", "hFDR.R"))



cLasso <- function(X, y, alpha = 0.05){
    lambda <- cali_lambda(X, y, alpha)$lambda
    
    fit <- glmnet::glmnet(X, y, lambda = lambda,
                          intercept = F, standardize = F,
                          standardize.response = F, family = "gaussian")
    selected <- which(fit$beta[, 1] != 0)
    
    return(list(selected = unname(selected)))
}


cali_lambda <- function(X, y, alpha, lambda_seq = NULL){
    n <- NROW(X)
    p <- NCOL(X)
    
    if(is.null(lambda_seq)){
        nlambda <- 100
        
        lambda_max <- max(abs(t(X) %*% y)) / n
        lambda_min <- min(lambda_max * 0.5, 1/n)
        k <- (0:(nlambda-1)) / nlambda
        lambda_seq <- lambda_max * (lambda_min/lambda_max)^k
    } else{
        nlambda <- length(lambda_seq)
    }
    
    hFDR <- hFDR.fixedX(X, y, tune_seq = lambda_seq,
                        filter_ind = round(nlambda*0.9), method = "lasso")
    
    lambda <- min(lambda_seq[hFDR <= alpha])
    
    return(list(lambda = lambda, hFDR = hFDR))
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

dir_trunc <- function(x, direction){
    if(length(x) == 0) return(NULL)
    
    # x <- sapply(x, function(e){
    #   if(direction * e > 0){e} else{direction * Inf}
    # })
    for(i in 1:length(x)){
        if(direction * x[i] <= 0){
            x[i] <- direction * Inf
        }
    }
    return(x)
}

solve_mat22 <- function(mat){
    a <- mat[1, 1]
    b <- mat[1, 2]
    c <- mat[2, 1]
    d <- mat[2, 2]
    det <- a*d-b*c
    inv <- matrix(c(d, -c, -b, a), ncol = 2) / det
    return(inv)
}

