library(glmnet)

cLasso <- function(X, y, alpha = 0.05){
    n <- NROW(X)
    p <- NCOL(X)
    
    mc_size <- 100
    nlambda <- 100
    
    lambda_max <- max(abs(t(X) %*% y)) / n
    lambda_min <- lambda_max / 2e3
    k <- (0:(nlambda-1)) / nlambda
    lambda_seq <- lambda_max * (lambda_min/lambda_max)^k
    
    DR <- sapply(1:p, function(j){
        y.samples <- y_condj_sample(y, X, j, mc_size)
        DP_j <- sapply(1:mc_size, function(mc_i){
            fit <- glmnet::glmnet(X, y.samples[, mc_i], lambda = lambda_seq,
                                  intercept = T, standardize = F,
                                  standardize.response = F, family = "gaussian")
            j_selected <- (fit$beta[j, ] != 0)
            selected_num <- colSums(fit$beta != 0)
            
            return(j_selected / pmax(selected_num, 1))
        })
        
        DR_j <- rowMeans(DP_j)
        
        return(DR_j)
    })
    DR <- rowSums(DR)
    
    lambda <- min(lambda_seq[DR <= alpha])
    
    fit <- glmnet::glmnet(X, y, lambda = lambda,
                          intercept = T, standardize = F,
                          standardize.response = F, family = "gaussian")
    selected <- which(fit$beta[, 1] != 0)
    
    return(list(selected = unname(selected)))
}



# naively sampling y conditional on Sj
y_condj_sample <- function(y, X, j, sample_size){
    
    n <- NROW(X)
    
    projj_y <- lm(formula = y ~ X[, -j] + 0)$fitted.values
    radius <- sqrt(sum(y^2) - sum(projj_y^2))
    
    y_sample <- matrix(rnorm(n * sample_size), nrow = n)
    
    y_sample <- y_sample - lm(formula = y_sample ~ X[, -j] + 0)$fitted.values
    y_sample <- scale(y_sample, center = FALSE, scale = sqrt(colSums(y_sample^2)) / radius)
    y_sample <- projj_y + y_sample
    
    return(y_sample)
}



# p <- 30; n <- 100; k <- 15
# X <- matrix(rnorm(n*p), n)
# nonzero <- sample(p, k)
# beta <- 1 * (1:p %in% nonzero)
# # beta <- 0.2 * (1:p %in% nonzero)
# y <- X %*% beta + rnorm(n) + 0
# print(which(1:p %in% nonzero))
# 
# # Basic usage
# result <- cLasso(X, y, alpha = 0.2)
# print(result$selected)
# 
# fdp <- function(selected) sum(beta[selected] == 0) / max(1, length(selected))
# fdp(result$selected)
# 
# tpp <- function(selected) sum(beta[selected] != 0) / max(1, sum(beta != 0))
# tpp(result$selected)