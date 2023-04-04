library(leaps)
library(glasso)

hFDR.fixedX <- function(X, y, tune_seq, filter_ind,
                        method = c("lasso", "FS"), mc_size = 50){
  mc_size <- 50
  
  method <- match.arg(method)
  p <- NCOL(X)
  
  if(method == "lasso"){
    lasso.fit <- glmnet::glmnet(X, y, lambda = tune_seq,
                                intercept = F, standardize = F,
                                standardize.response = F, family = "gaussian")
    hFDR <- lasso.hFDR.simple(X, y, lasso.fit)$hFDR
  } else if(method == "FS"){
    tune_filter <- max(tune_seq[filter_ind], floor(p*0.8))
    filter.res <- select_variables(X, y, tune_filter, method)
    hNulls <- which(!filter.res)
    
    DRj_tune <- matrix(0, nrow = p, ncol = length(tune_seq))
    filter_prob_unselect <- rep(1, p)
    
    for(j in hNulls){
      DPj_tune <- matrix(0, nrow = length(tune_seq), ncol = mc_size)
      y.samples <- sampler.fixedX(X, y, j, mc_size)
      
      for(mc_i in 1:mc_size){
        X.sample <- X
        y.sample <- y.samples[, mc_i]
        res.sample <- select_variables(X.sample, y.sample, tune_seq, method)
        DPj_tune[, mc_i] <- res.sample[j, ] / pmax(1, colSums(res.sample))
      }
      
      DRj_tune[j, ] <- rowMeans(DPj_tune)
      filter_prob_unselect[j] <- sum(DPj_tune[filter_ind, ] == 0) / mc_size
    }
    hFDR <- colSums(DRj_tune / filter_prob_unselect)
  }
  
  return(hFDR)
}


hFDR.modelX <- function(X, y, Xcov.true, tune_seq, filter_ind,
                        method = c("lasso", "FS"), mc_size = 50){
  
  method <- match.arg(method)
  p <- NCOL(X)
  
  tune_filter <- if(method == "FS"){
    max(tune_seq[filter_ind], floor(p*0.8))
  } else{
    tune_seq[filter_ind]
  }
  filter.res <- select_variables(X, y, tune_filter, method)
  hNulls <- which(!filter.res)
  
  DRj_tune <- matrix(0, nrow = p, ncol = length(tune_seq))
  filter_prob_unselect <- rep(1, p)
  
  for(j in hNulls){
    DPj_tune <- matrix(0, nrow = length(tune_seq), ncol = mc_size)
    X.samples <- sampler.modelX(X, Xcov.true, j, mc_size)
    
    for(mc_i in 1:mc_size){
      X.sample <- X
      X.sample[, j] <- X.samples[, mc_i]
      y.sample <- y
      res.sample <- select_variables(X.sample, y.sample, tune_seq, method)
      DPj_tune[, mc_i] <- res.sample[j, ] / pmax(1, colSums(res.sample))
    }
    
    DRj_tune[j, ] <- rowMeans(DPj_tune)
    filter_prob_unselect[j] <- sum(DPj_tune[filter_ind, ] == 0) / mc_size
  }
  hFDR <- colSums(DRj_tune / filter_prob_unselect)
  
  return(hFDR)
}


hFDR.glasso <- function(X, tune_seq, filter_ind,
                        method = c("glasso"), mc_size = 50){
  
  method <- match.arg(method)
  p <- NCOL(X)
  pair_num <- p*(p-1)/2
  
  filter.res <- select_variables(X, y = NA, tune_seq[filter_ind], method)
  hNulls <- matrix(F, p, p)
  hNulls[upper.tri(hNulls)] <- !filter.res
  
  DRj_tune <- matrix(0, nrow = pair_num, ncol = length(tune_seq))
  filter_prob_unselect <- rep(1, pair_num)
  
  for(j in 2:p){
    for(i in 1:(j-1)){
      if(!hNulls[i,j]) next
      
      DPj_tune <- matrix(0, nrow = length(tune_seq), ncol = mc_size)
      X.samples <- sampler.glassoX(X, i, j, mc_size)
      
      for(mc_i in 1:mc_size){
        X.sample <- X
        X.sample[, i] <- X.samples$Xi[, mc_i]
        X.sample[, j] <- X.samples$Xj[, mc_i]
        res.sample <- select_variables(X.sample, NA, tune_seq, method)
        DPj_tune[, mc_i] <- res.sample[(j-1)*(j-2)/2 + i, ] / pmax(1, colSums(res.sample))
      }
      
      DRj_tune[(j-1)*(j-2)/2 + i, ] <- rowMeans(DPj_tune)
      filter_prob_unselect[(j-1)*(j-2)/2 + i] <- sum(DPj_tune[filter_ind, ] == 0) / mc_size
    }
  }
  hFDR <- colSums(DRj_tune / filter_prob_unselect)
  
  return(hFDR)
}

select_variables <- function(X, y, tune_seq, method = c("lasso", "FS", "glasso")){
  method <- match.arg(method)
  
  if(method == "lasso"){
    res <- glmnet::glmnet(X, y, lambda = tune_seq,
                          intercept = F, standardize = F,
                          standardize.response = F, family = "gaussian")
    res <- (res$beta != 0)
  } else if(method == "FS"){
    res <- regsubsets(x = X, y = y, nvmax = max(tune_seq), method = "forward")
    tryCatch({
      res <- t(unname(summary(res)$which[tune_seq, -1]))
    }, error = function(e){browser()})
  } else if(method == "glasso"){
    S <- var(X)
    res <- sapply(tune_seq, function(rho){
      inv.Sigma.est <- glasso(S, rho = rho)$wi
      inv.Sigma.est[upper.tri(inv.Sigma.est)] != 0
    })
  } else{
    stop()
  }
  
  return(res)
}


sampler.fixedX <- function(X, y, j, sample_size){
  sampler.core(y, X[, -j], sample_size)
}

sampler.modelX <- function(X, Sigma, j, sample_size){
  n <- NROW(X)
  
  trans.mat <- Sigma[j, -j] %*% solve(Sigma[-j, -j])
  mean.cond <- c(trans.mat %*% t(X[, -j]))
  cov.cond <- c(Sigma[j, j] - trans.mat %*% Sigma[-j, j])
  
  Xj.sample <- mean.cond + sqrt(cov.cond) * matrix(rnorm(n * sample_size), nrow = n)
}

sampler.glassoX <- function(X, i, j, sample_size){
  Xi.samples <- sampler.core(X[, i], X[, -c(i, j)], sample_size)
  Xj.samples <- sampler.core(X[, j], X[, -c(i, j)], sample_size)
  
  return(list(Xi = Xi.samples, Xj = Xj.samples))
}

sampler.core <- function(var, basis, sample_size){
  n <- length(var)
  
  var.proj <- lm(var ~ basis + 0)$fitted.values
  
  radius <- sqrt(sum(var^2) - sum(var.proj^2))
  
  var.sample <- matrix(rnorm(n * sample_size), nrow = n)
  
  var.sample.proj <-lm(var.sample ~ basis + 0)$fitted.values
  
  var.sample <- var.sample - var.sample.proj
  var.sample <- scale(var.sample, center = FALSE, scale = sqrt(colSums(var.sample^2)) / radius)
  var.sample <- var.proj + var.sample
  
  return(var.sample)
}