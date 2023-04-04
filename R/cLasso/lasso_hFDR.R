library(glmnet)

lasso.hFDR <- function(X, y, lasso.cv){
  if(!all(class(lasso.cv) == "cv.glmnet")){
    stop()
  }
  
  n <- NROW(X)
  p <- NCOL(X)
  data.pack <- process_data(X, y)
  
  sigma_hat <- sqrt(data.pack$RSS_X / (n-p))
  
  nlambda <- length(lasso.cv$lambda)
  nsample <- 20
  lambda0.index <- which.min(abs(lasso.cv$lambda - sigma_hat))
  lambda.end <- max(lambda0.index, lasso.cv$index + nlambda / nsample)
  lambda.end <- min(nlambda, round(lambda.end))
  
  ind.sample <- round(seq(from = 1, to = lambda.end, length.out = nsample))
  ind.sample <- sort(unique(c(ind.sample, lasso.cv$index)))
  lambda <- lasso.cv$lambda[ind.sample]
  lasso.fit <- list(lambda = lambda, beta = lasso.cv$glmnet.fit$beta[, ind.sample])
  
  glmnet.pack <- pack_glmnet(X, y, lasso.fit)
  data.pack <- process_data(X, y)
  
  hFDR.res <- calc_lasso_hFDR(glmnet.pack, data.pack, couple = T,
                              denoise = F, var_est = T, mc_size = 30,
                              beta_hat.lasso = lasso.cv$glmnet.fit$beta[, lasso.cv$index[2]],
                              lambda_star = lasso.cv$lambda.1se * n)
  
  hFDR.res <- structure(list(call = match.call(),
                             X = X,
                             y = y,
                             glmnet.cv = lasso.cv,
                             lambda = lambda,
                             hFDR = hFDR.res$hFDR_lambda,
                             hFDR.std = hFDR.res$deviate.est,
                             hFDR.low = hFDR.res$hFDR_lambda - hFDR.res$deviate.est,
                             hFDR.up = hFDR.res$hFDR_lambda + hFDR.res$deviate.est,
                             name = "estimated FDR"),
                        class = 'lasso.hFDR')
  
  return(hFDR.res)
}

lasso.hFDR.simple <- function(X, y, lasso.fit){
  glmnet.pack <- pack_glmnet(X, y, lasso.fit)
  data.pack <- process_data(X, y)
  
  hFDR.res <- calc_lasso_hFDR(glmnet.pack, data.pack, couple = T,
                              denoise = F, var_est = F, mc_size = 30,
                              beta_hat.lasso = NA,
                              lambda_star = NA)
  
  hFDR.res <- structure(list(lambda = lasso.fit$lambda,
                             hFDR = hFDR.res$hFDR_lambda,
                             hFDRj = hFDR.res$DRj_lambda),
                        class = 'lasso.hFDR.simple')
  
  return(hFDR.res)
}

lasso.FDR.select <- function(X, y, lasso.fit, alpha){
  lambdas <- lasso.fit$lambda
  
  hFDR.obj <- lasso.hFDR.simple(X, y, lasso.fit)
  
  lambda.hat.index <- max(which(hFDR.obj$hFDR <= alpha))
  if(abs(lambda.hat.index) == Inf) stop("lambda range too narrow")
  
  lambda.cand.index <- max(which(hFDR.obj$hFDR <= 1.5 * alpha))
  candidates <- which(lasso.fit$beta[, lambda.cand.index] != 0)
  
  first_nonzero <- function(x) match(T, abs(x) > 0) # NA if all(x==0)
  time_of_entry <- apply(lasso.fit$beta, 1, first_nonzero)
  names(time_of_entry) <- NULL
  time_of_entry[is.na(time_of_entry)] <- nlambda + 1
  
  mc_size <- 50
  print(candidates)
  preselected <- sapply(candidates, function(j){
    print(j)
    y_samples <- y_condj_sample(y, X, j, mc_size)
    
    bj.mc <- rep(NA, mc_size)
    DPj.mc <- rep(NA, mc_size)
    for(mc_i in 1:mc_size){
      lasso.fit.mc <- glmnet::glmnet(X, y_samples[, mc_i], lambda = lambdas,
                                     intercept = F, standardize = F,
                                     standardize.response = F, family = "gaussian")
      hFDR.obj.mc <- lasso.hFDR.simple(X, y_samples[, mc_i], lasso.fit.mc)
      
      lambda.hat.mc.index <- max(which(hFDR.obj.mc$hFDR <= alpha))
      if(abs(lambda.hat.mc.index) == Inf) stop("lambda range too narrow")
      
      tvals <- lm_to_t(y_samples[, mc_i], X)
      pvals <- pvals_t(tvals$tvals, n-p, side = "two")
      
      bj.mc[mc_i] <- (pvals[j] > 0.5) *2 * hFDR.obj$hFDRj[j, lambda.hat.mc.index]
      
      time_of_entry.mc <- apply(lasso.fit.mc$beta, 1, first_nonzero)
      names(time_of_entry.mc) <- NULL
      time_of_entry.mc[is.na(time_of_entry.mc)] <- nlambda + 1
      
      DPj.mc[mc_i] <- (time_of_entry.mc[j] <= time_of_entry[j]) / 
        length(union(which(lasso.fit.mc$beta[, lambda.hat.mc.index] != 0), j))
    }
    # browser()
    
    return(mean(DPj.mc) <= mean(bj.mc))
  })
  preselected <- candidates[preselected]
  
  R.hat <- which(lasso.fit$beta[, lambda.hat.index] != 0)
  Rj.hat <- sapply(preselected, function(j){
    length(union(R.hat, j))
  })
  rand.prune <- (length(preselected) > 0) && any(length(preselected) < Rj.hat)
  if(rand.prune){
    u <- runif(length(preselected))
    R <- max(which(sapply(1:length(preselected), function(r){
      sum(u <= r / Rj.hat) >= r
    })))
    selected <- preselected[u <= R / Rj.hat]
  } else{
    selected <- preselected
  }
  
  return(list(naive.selected = R.hat, pre.selected = preselected, 
              rand.prune = rand.prune, selected = selected))
}

pack_glmnet <- function(X, y, lasso.fit, tol = 1e-7){
  n <- NROW(X)
  p <- NCOL(X)
  
  Sigma <- t(X) %*% X
  Xy <- c(t(y) %*% X)
  
  lambdas <- lasso.fit$lambda * n
  nlambda <- length(lambdas)
  
  selected.list <- list()
  Sigma.selected.inv.list <- list()
  pre_selected <- NA
  for(lambda_i in 1:nlambda){
    selected <- sort(which(abs(lasso.fit$beta[, lambda_i]) > tol))
    
    if(!identical(pre_selected, selected)){
      Sigma.selected.inv <- if(length(selected) > 0){
        list(solve(Sigma[selected, selected]))
      } else{ list(NULL) }
    }
    
    Sigma.selected.inv.list[lambda_i] <- Sigma.selected.inv
    selected.list[[lambda_i]] <- selected
    
    pre_selected <- selected
  }
  
  return(list(Sigma = Sigma, Xy = Xy,
              lambdas = lambdas, betas = lasso.fit$beta,
              selected.list = selected.list,
              Sigma.selected.inv.list = Sigma.selected.inv.list,
              tol = tol))
}

calc_lasso_hFDR <- function(glmnet.pack, data.pack, couple = T,
                            denoise = F, var_est = F, mc_size = 10,
                            beta_hat.lasso = NA, lambda_star = NA){
  X <- data.pack$X
  y <- data.pack$y
  tol <- glmnet.pack$tol
  
  n <- NROW(X)
  p <- NCOL(X)
  
  sigma_hat <- sqrt(data.pack$RSS_X / (n-p))
  
  lambdas <- glmnet.pack$lambdas
  nlambda <- length(lambdas)
  
  DRj_lambda <- matrix(NA, nrow = p, ncol = nlambda)
  DRj_lambda.org <- matrix(NA, nrow = p, ncol = nlambda)
  
  lambda0.index <- which.min(abs(lambdas - sigma_hat))
  lambda0_selected <- glmnet.pack$selected.list[[lambda0.index]]
  lambda0_prob_unselect <- rep(1, p)
                                       
  
  tvals <- data.pack$vjy_obs / sqrt(data.pack$RSS_X / (n-p))
  pvals <- pvals_t(tvals, n-p, side = "two")

  # # caution ############
  # lambda0_selected <- which(pvals <= 0.5)

  if(denoise){
    n_denoise <- ceiling(p/2)
    # null_hypo <- which(pvals > 0.5)
    # if(length(null_hypo) >= n_denoise){
    #   denoise_expc_set <- sample(null_hypo, n_denoise)
    # } else{
    #   denoise_expc_set <- c(null_hypo, sample(setdiff(1:p, null_hypo), n_denoise - length(null_hypo)))
    # }
    denoise_expc_set <- sample(1:p, n_denoise)
    denoise_expc_set <- sort(denoise_expc_set)
    denoise_expcs <- matrix(0, nrow = p, ncol = nlambda)
  }
    
  for(j in 1:p){
    
    if(denoise && j %in% denoise_expc_set){
      y.samples <- y_condj_sample(y, X, j, mc_size, data.pack)
      mc.pack <- list()
      
      X.pack <- list(X = data.pack$X, Q_X = data.pack$Q_X, vj_mat = data.pack$vj_mat,
                     level_score = data.pack$level_score)
      for(mc_i in 1:mc_size){
        data.pack.mc <- process_y(X.pack, y.samples[, mc_i])
        mc.pack[[mc_i]] <- list(data.pack = data.pack.mc,
                                beta_lasso = matrix(NA, nrow = p, ncol = nlambda),
                                selected.list = list(),
                                Sigma.selected.inv.list = list())
      }
    }
    
    for(lambda_i in 1:nlambda){
      selected <- glmnet.pack$selected.list[[lambda_i]]
      lasso.homopath <- NULL
      
      lasso.pack <- list(Xy = glmnet.pack$Xy,
                         beta_lasso = glmnet.pack$betas[, lambda_i],
                         selected = selected,
                         Sigma = glmnet.pack$Sigma,
                         Sigma.selected.inv = glmnet.pack$Sigma.selected.inv.list[[lambda_i]],
                         lambda = lambdas[lambda_i],
                         tol = tol)
      # tryCatch({
      #   lasso.pack <- list(Xy = glmnet.pack$Xy,
      #                      beta_lasso = glmnet.pack$betas[, lambda_i],
      #                      selected = selected,
      #                      Sigma = glmnet.pack$Sigma,
      #                      Sigma.selected.inv = glmnet.pack$Sigma.selected.inv.list[[lambda_i]],
      #                      lambda = lambdas[lambda_i],
      #                      tol = tol)
      # }, error = function(e){browser()})
      Xjy_range <- c(data.pack$Xy_bound[, j])
      
      # if(j %in% lambda0_selected){
      if(F){
        DRj_lambda[j, lambda_i] <- 0
      }
      else{
        if(couple){
          lasso.homopath <- lasso_homotopy(lasso.pack, j, Xjy_range)
          res <- calc_DRj(lasso.homopath, j, Xjy_range, data.pack, tol = tol)
          
          DRj_lambda[j, lambda_i] <- res$DRj
          prob_j_unselect <- 1 - res$Evj
        } else{
          prob_j_unselect <- calc_prob_unselect(lasso.pack, j, data.pack)
          DRj_lambda[j, lambda_i] <- (1 - prob_j_unselect) / max(length(selected), 1)
        }
        
        if(lambda_i == lambda0.index){
          lambda0_prob_unselect[j] <- prob_j_unselect
        }
        
        # # caution ############
        # DRj_lambda.org[j, lambda_i] <- DRj_lambda[j, lambda_i]
        # if(j %in% lambda0_selected){
        #   DRj_lambda[j, lambda_i] <- 0
        # }
      }
      
      if(denoise && j %in% denoise_expc_set){
        if(is.null(lasso.homopath)){
          lasso.homopath <- lasso_homotopy(lasso.pack, j, Xjy_range)
        }
        
        for(mc_i in 1:mc_size){
          Xjy.sample <- sum(X[, j] * y.samples[, mc_i])
          lasso.sample <- homotopy_fit(Xjy.sample, lasso.homopath)
          
          mc.pack[[mc_i]]$beta_lasso[, lambda_i] <- lasso.sample$beta_lasso
          mc.pack[[mc_i]]$selected.list[[lambda_i]] <- lasso.sample$selected
          mc.pack[[mc_i]]$Sigma.selected.inv.list[[lambda_i]] <- lasso.sample$Sigma.selected.inv
        }
      }
    }
    
    if(denoise && j %in% denoise_expc_set){
      denoise_expcs[j, ] <- rowMeans(sapply(1:mc_size, function(mc_i){
        glmnet.pack.mc <- list(Sigma = glmnet.pack$Sigma,
                               Xy = c(t(y.samples[, mc_i]) %*% X),
                               lambdas = lambdas, betas = mc.pack[[mc_i]]$beta_lasso,
                               selected.list = mc.pack[[mc_i]]$selected.list,
                               Sigma.selected.inv.list = mc.pack[[mc_i]]$Sigma.selected.inv.list,
                               tol = tol)
        data.pack.mc <- mc.pack[[mc_i]]$data.pack
        denoiser.mc <- calc_lasso_hFDR(glmnet.pack.mc, data.pack.mc, couple = F,
                                       denoise = F)$hFDR_lambda
        return(denoiser.mc)
      }))
    }
    
  }
  
  # # caution ############
  # lambda0_prob_unselect <- rep(0.5, p)
  
  DRj_lambda <- DRj_lambda / lambda0_prob_unselect
  
  hFDR_lambda <- colSums(DRj_lambda)

  if(denoise){
    denoiser <- calc_lasso_hFDR(glmnet.pack, data.pack, couple = F,
                                denoise = F)$hFDR_lambda
    denoise_expcs <- colSums(denoise_expcs) / n_denoise

    hFDR_lambda <- hFDR_lambda - denoiser + denoise_expcs
  }
  
  # if(denoise){
  #   denoiser <- calc_lasso_hFDR(glmnet.pack, data.pack, couple = F,
  #                               denoise = F)$hFDR_lambda
  #   corrector <- matrix(0, nrow = p, ncol = nlambda)
  #   for(j in denoise_expc_set){
  #     corrector[j, ] <- denoiser
  #   }
  #   corrector <- (corrector - denoise_expcs) / n_denoise
  #   
  #   DRj_lambda <- pmax(DRj_lambda - corrector, 0)
  # }
  # hFDR_lambda <- colSums(DRj_lambda)
  
  
  if(var_est){
    beta_hat.ols <- data.pack$vjy_obs * sqrt(data.pack$level_score)
    # beta_hat <- beta_hat.ols
    
    # lambda.index <- 8
    # lambda_star <- lambdas[lambda.index]
    # beta_hat.lasso <- c(glmnet.pack$betas[, lambda.index])
    
    # beta_hat <- beta_hat.lasso
    # beta_hat[beta_hat != 0] <- beta_hat[beta_hat != 0] + sign(beta_hat[beta_hat != 0]) * lambda_star
    
    # beta_hat <- rep(0, p)
    # beta_hat[beta_hat.lasso != 0] <- beta_hat.ols[beta_hat.lasso != 0]
    
    beta_hat <- rep(0, p)
    if(any(beta_hat.lasso != 0)){
      beta_hat[beta_hat.lasso != 0] <- lm(y ~ X[, beta_hat.lasso != 0] + 0)$coefficients
    }
    
    X.pack <- list(X = data.pack$X, Q_X = data.pack$Q_X, vj_mat = data.pack$vj_mat,
                   level_score = data.pack$level_score)
    
    hFDR.samples <- matrix(NA, nlambda, mc_size)
    FDP.samples <- matrix(NA, nlambda, mc_size)
    H0.boot <- beta_hat == 0
    set.seed(-100)
    
    for(boot_i in 1:mc_size){
    # for(boot_i in 1:n){
      X.boot <- X
      y.boot <- X %*% beta_hat + rnorm(n, sd = sigma_hat)
      data.pack.boot <- process_y(X.pack, y.boot)
      
      # bs.samples <- sample(1:n, n, replace = T)
      # # bs.samples <- (1:n)[-boot_i]
      # X.boot <- X[bs.samples, ]
      # y.boot <- y[bs.samples]
      # data.pack.boot <- process_data(X.boot, y.boot)
      
      lasso.fit.boot <- glmnet::glmnet(X.boot, y.boot, lambda = lambdas / n,
                                       intercept = F, standardize = F,
                                       standardize.response = F, family = "gaussian")
      glmnet.pack.boot <- pack_glmnet(X.boot, y.boot, lasso.fit.boot)
      
      hFDR.boot <- calc_lasso_hFDR(glmnet.pack.boot, data.pack.boot, couple = F,
                                   denoise = denoise, var_est = F, mc_size = mc_size)$hFDR_lambda
      FDP.boot <- calc_lasso_FDP(lasso.fit.boot, H0.boot)
      
      hFDR.samples[, boot_i] <- hFDR.boot
      FDP.samples[, boot_i] <- FDP.boot
    }
    deviate.est <- calc_deviate(hFDR.samples, FDP.samples)
    # deviate.est <- calc_deviate(hFDR.samples, FDP.samples) * abs(qt(0.16, df = n-1))
  } else{
    deviate.est <- rep(NA, nlambda)
    beta_hat.ols <- rep(NA, p)
    beta_hat.lasso <- rep(NA, p)
    beta_hat <- rep(NA, p)
    sigma_hat <- NA
  }
  
  return(list(hFDR_lambda = hFDR_lambda, DRj_lambda = DRj_lambda.org,
              deviate.est = deviate.est,
              beta_hat.ols = beta_hat.ols, beta_hat.lasso = beta_hat.lasso,
              beta_hat = beta_hat, sigma_hat = sigma_hat, pvals = pvals))
}

calc_DRj <- function(lasso.homopath, j, Xjy_range, data.pack, tol = 1e-7){
  res_norm2 <- data.pack$RSS_Xnoj[j]
  df <- data.pack$n - data.pack$p
  trans <- data.pack$trans
  
  n_nodes <- length(lasso.homopath$Xjy_nodes)
  mid_beta <- (lasso.homopath$beta_at_nodes[, 1:(n_nodes-1)] + 
    lasso.homopath$beta_at_nodes[, 2:n_nodes]) / 2
  mid_beta <- matrix(mid_beta, ncol = n_nodes-1)
  DPj <- (abs(mid_beta[j, ]) > tol) / pmax(colSums(abs(mid_beta) > tol), 1)
  
  trunc.low <- sum(lasso.homopath$Xjy_nodes <= min(Xjy_range))
  trunc.up <- sum(lasso.homopath$Xjy_nodes >= max(Xjy_range))
  if(trunc.low < 1 || trunc.up < 1) browser()
  main_trunk <- if(trunc.low+1 <= n_nodes-trunc.up){
    lasso.homopath$Xjy_nodes[(trunc.low+1):(n_nodes-trunc.up)]
  } else { NULL }
  # int_nodes <- c(min(Xjy_range), main_trunk, max(Xjy_range))
  int_nodes <- main_trunk
  int_nodes <- Xjy_to_vjy(int_nodes, trans, j)
  DPj <- DPj[trunc.low:(n_nodes-trunc.up)]
  
  CDF <- c(trans$Xv[j] < 0, vjy_CDF(int_nodes, res_norm2, df), trans$Xv[j] > 0)
  DRj <- sum(abs(diff(CDF)) * DPj)
  Evj <- sum(abs(diff(CDF)) * (DPj > 0))
  
  return(list(DRj = DRj, Evj = Evj))
}

calc_prob_unselect <- function(lasso.pack, j, data.pack){
  p <- data.pack$p
  n <- data.pack$n
  df <- n - p
  
  Xjy_range <- c(data.pack$Xy_bound[, j])
  
  if(!(j %in% lasso.pack$selected)){
    Xjy_nodes <- c()
    for(direction in c(1, -1)){
      res_at_node <- lasso_homotopy_step(lasso.pack, j, Xjy_range, direction)
      Xjy_nodes <- c(Xjy_nodes, res_at_node$Xjy)
    }
    Xjy_nodes[Xjy_nodes <= min(Xjy_range)] <- min(Xjy_range)
    Xjy_nodes[Xjy_nodes >= max(Xjy_range)] <- max(Xjy_range)
    int_nodes <- sort(Xjy_to_vjy(Xjy_nodes, data.pack$trans, j))
    prob_unselect <- abs(diff(vjy_CDF(int_nodes, data.pack$RSS_Xnoj[j], df)))
  } else{
    lasso.homopath <- lasso_homotopy(lasso.pack, j, Xjy_range)
    Evj <- calc_DRj(lasso.homopath, j, Xjy_range, data.pack,
                    tol = lasso.pack$tol)$Evj
    prob_unselect <- 1 - Evj
  }
  
  return(prob_unselect)
}

plot.lasso.hFDR <- function(hFDR.obj, H0, sign.lambda = -1,...){
  xlab <- expression(Log(lambda))
  if(sign.lambda < 0) xlab <- paste("-", xlab, sep = "")
  
  nlambda <- length(hFDR.obj$glmnet.cv$lambda)
  lasso.FDP <- sapply(1:nlambda, function(lambda_index){
    lasso.selected <- which(lasso.cv$glmnet.fit$beta[, lambda_index] != 0)
    false_discover <- length(intersect(lasso.selected, which(H0)))
    false_discover / max(length(lasso.selected), 1)
  })
  
  
  plot.range <- range(hFDR.obj$hFDR.low, hFDR.obj$hFDR.up)
  plot.args <- list(x = sign.lambda * log(hFDR.obj$lambda),
                    y = hFDR.obj$hFDR,
                    xlim = range(sign.lambda*log(hFDR.obj$glmnet.cv$lambda)),
                    ylim = plot.range,
                    xlab = xlab, ylab = "estimated FDR and scaled CV MSE", type = "n")
  new.args <- list(...)
  if(length(new.args)) plot.args[names(new.args)] <- new.args
  
  do.call("plot", plot.args)
  lines(x = sign.lambda * log(hFDR.obj$glmnet.cv$lambda),
        y = (hFDR.obj$glmnet.cv$cvm - min(hFDR.obj$glmnet.cv$cvm)) / abs(diff(range(hFDR.obj$glmnet.cv$cvm))) * abs(diff(plot.range)) + min(plot.range),
        col = "red")
  lines(x = sign.lambda * log(hFDR.obj$glmnet.cv$lambda),
        y = lasso.FDP,
        col = "black")
  error.bars(sign.lambda * log(hFDR.obj$lambda),
             hFDR.obj$hFDR.up, hFDR.obj$hFDR.low,
             width = 0.01, col = "darkgrey")
  points(sign.lambda*log(hFDR.obj$lambda), hFDR.obj$hFDR,
         pch = 20, col = "blue")
  axis(side = 3, at = sign.lambda*log(hFDR.obj$glmnet.cv$lambda),
       labels = paste(hFDR.obj$glmnet.cv$nz), tick = FALSE, line = 0)
  abline(v = sign.lambda * log(hFDR.obj$glmnet.cv$lambda.min), lty = 3)
  abline(v = sign.lambda * log(hFDR.obj$glmnet.cv$lambda.1se), lty = 3)
  invisible()
}

plot_path <- function(hFDR.obj, sign.lambda = -1,...){
  xlab <- expression(Log(lambda))
  if(sign.lambda < 0) xlab <- paste("-", xlab, sep = "")
  
  nbeta <- NCOL(hFDR.obj$glmnet.cv$glmnet.fit$beta)
  plot.args <- list(x = sign.lambda * log(hFDR.obj$lambda),
                    xlim = range(sign.lambda*log(hFDR.obj$glmnet.cv$lambda)),
                    ylim = range(as.vector(hFDR.obj$glmnet.cv$glmnet.fit$beta)),
                    xlab = xlab, ylab = "coefficients", type = "n")
  new.args <- list(...)
  if(length(new.args)) plot.args[names(new.args)] <- new.args
  
  do.call("plot", plot.args)
  for(beta_i in 1:nbeta){
    lines(x = sign.lambda * log(hFDR.obj$glmnet.cv$lambda),
          y = as.vector(hFDR.obj$glmnet.cv$glmnet.fit$beta[beta_i, ]),
          col = "grey")
  }
  invisible()
}



error.bars <- function(x, upper, lower, width = 0.02, ...) {
  xlim <- range(x)
  barw <- diff(xlim) * width
  segments(x, upper, x, lower, ...)
  segments(x - barw, upper, x + barw, upper, ...)
  segments(x - barw, lower, x + barw, lower, ...)
  range(upper, lower)
}











process_X <- function(X){
  n <- NROW(X)
  p <- NCOL(X)
  
  QR <- qr(X)
  
  pivot_back <- sort.list(QR$pivot)
  
  Q_X <- qr.Q(QR, complete = F)
  R_X <- qr.R(QR, complete = F)[, pivot_back]
  
  # compute the a basic matrix needed by cknockoff
  # compute vj = unit(X_{j.-j}), X_{j.-j} = X_j orthogonal projected onto X_{-j}. time O(p^3)
  vj_mat <- sapply(1:p, function(j){
    # Q_X[, j] = X_j orthogonal projected onto X_{1:j-1}
    # X_{j.-j} = Q_X[, j] orthogonal projected onto S, S:=(X_{j+1:p} orthogonal projected onto X_{1:j-1})
    #          <=> find a vector in span(X_{j:p}) that is perpendicular to S
    # "coord" is the coordinate of such a vector under the basis Q_X[, j:p]
    coord <- forwardsolve(t(R_X[j:p,j:p]), c(1, rep(0, p-j)))
    vj <- Q_X[, j:p] %*% matrix(coord, nrow = p-j+1)
    vj <- vj / sqrt(sum(vj^2))
  })
  
  return(list(X = X, Q_X = Q_X, vj_mat = vj_mat, level_score = diag(solve(t(X)%*%X))))
}

process_y <- function(X.pack, y, ep = 1e-4){
  X <- X.pack$X
  Q_X <- X.pack$Q_X
  vj_mat <- X.pack$vj_mat
  
  n <- NROW(X)
  p <- NCOL(X)
  
  vjy_obs <- c(matrix(y, nrow=1) %*% vj_mat)
  RSS_X <- abs(sum(y^2) - sum((matrix(y, nrow=1) %*% Q_X)^2))
  RSS_Xnoj <- RSS_X + vjy_obs^2
  
  Xv <- colSums(X * vj_mat)
  X_perpv_y <- c(t(y) %*% X) - Xv * vjy_obs
  
  trans <- list(Xv = Xv, X_perpv_y = X_perpv_y)
  
  vy_bound <- pmax(abs(vjy_quantile(ep, RSS_Xnoj, df = n-p)), abs(vjy_obs))
  Xy_bound <- sapply(1:p, function(j){
    sort(c(vjy_to_Xjy(-vy_bound[j], trans, j),
           vjy_to_Xjy(vy_bound[j], trans, j)))
  })
  
  return(list(X = X, y = y, Q_X = Q_X, vj_mat = vj_mat, level_score = X.pack$level_score,
              RSS_X = RSS_X, RSS_Xnoj = RSS_Xnoj, trans = trans,
              Xy_bound = Xy_bound, vjy_obs = vjy_obs,
              n = n, p = p))
}

process_data <- function(X, y, ep = 1e-4){
  X.pack <- process_X(X)
  data.pack <- process_y(X.pack, y, ep)
  
  return(data.pack)
}

Xjy_to_vjy <- function(Xjy, trans, j){
  (Xjy - trans$X_perpv_y[j]) / trans$Xv[j]
}

vjy_to_Xjy <- function(vjy, trans, j){
  trans$Xv[j] * vjy + trans$X_perpv_y[j]
}

# the distribution function of Xj^perp y conditional on Sj. See 11.summary Thm 3.1
vjy_CDF <- function(vjy, res_norm2, df){
  return(pt(vjy * sqrt(df) / sqrt(pmax(res_norm2 - vjy^2, 0)), df = df))
}

vjy_quantile <- function(prob, res_norm2, df){
  t_quantile <- qt(prob, df)
  if(abs(t_quantile) < Inf){
    vjy <- sign(t_quantile) * sqrt(t_quantile^2 * res_norm2 / (t_quantile^2 + df))
  } else{
    vjy <- sign(t_quantile) * sqrt(res_norm2)
  }
  return(vjy)
}







