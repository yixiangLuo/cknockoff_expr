library(here)

library(tidyverse)
library(knockoff)
library(cknockoff)

library(foreach)
library(doParallel)

source(here("R", "methods.R"))
source(here("R", "utils.R"))
source(here("R", "plot.R"))


experiment <- "scale_n_expr"

X_types <- c("IID_Normal", "MCC")
posit_type <- "random"

alpha <- 0.05
p <- 100
n_vec <- p * c(3, 5, 10, 20)
pi1 <- 10 / p
n_rounds <- rep(50, length(n_vec))  # 7 physical cores used

target <- 0.5
target_at_alpha <- 0.2

n_cores <- 7

X_seed <- 2021
noise <- quote(rnorm(n))

knockoffs.name <- "ckn.create.fixed"
statistic.name <- "stat.glmnet_coefdiff_tiebreak"

knockoffs <- get(knockoffs.name)
statistic <- get(statistic.name)

runtime.result <- lapply(X_types, function(X_type){
  print(X_type)
  
  runtime_vec <- sapply(1:length(n_vec), function(iter){
    n <- n_vec[iter]
    
    print(n)
    
    n_round <- n_rounds[iter]
    
    random_X.data <- list(random_X = T,
                          n = n, p = p, X_type = X_type, sample_num = 5)
    mu1 <- BH_lm_calib(X = NA, random_X.data, pi1, noise,
                       posit_type, 1, side = "two", nreps = 50,
                       alpha = target_at_alpha, target = target, n_cores = 14)
    
    beta <- genmu(p, pi1, mu1, posit_type, 1)
    
    H0 <- beta == 0
    
    print(paste0(c("alt:", which(!H0)), collapse = " "))
    
    registerDoParallel(n_cores)
    
    runtime <- foreach(iter = 1:n_round) %dopar% {
    # runtime <- sapply(1:n_round, function(iter){
      
      X <- gene_X(X_type, n, p, iter)$X
      
      y <- X %*% beta + eval(noise)
      
      time_start <- Sys.time()
      kn.result <- knockoff.filter(X, y, knockoffs = knockoffs,
                                   statistic = statistic, fdr = alpha)
      time_end <- Sys.time()
      
      kn.runtime <- difftime(time_end, time_start, units = "secs")[[1]]
      
      time_start <- Sys.time()
      ckn.result <- cknockoff(X, y,
                              intercept = F,
                              statistic = statistic,
                              alpha = alpha,
                              n_cores = 1,
                              Rstar_refine = F)
      time_end <- Sys.time()
      
      ckn.runtime <- difftime(time_end, time_start, units = "secs")[[1]]
      
      time_start <- Sys.time()
      ckn_star.result <- cknockoff(X, y,
                                   intercept = F,
                                   statistic = statistic,
                                   alpha = alpha,
                                   n_cores = 1,
                                   Rstar_refine = T)
      time_end <- Sys.time()
      
      ckn_star.runtime <- difftime(time_end, time_start, units = "secs")[[1]]
      
      print(paste0(c("kn:", kn.result$selected), collapse = " "))
      print(paste0(c("ckn:", ckn.result$selected), collapse = " "))
      print(paste0(c("ckn*:", ckn_star.result$selected), collapse = " "))
      
      runtime_sample <- c(kn.runtime, ckn.runtime, ckn_star.runtime)
      names(runtime_sample) <- c("knockoff", "cKnockoff", "cKnockoff_STAR")
      
      return(runtime_sample)
    }
    runtime <- do.call(cbind, runtime)
    runtime_avg <- rowMeans(runtime)
    
    return(runtime_avg)
  })
  
  runtime <- as.data.frame(t(runtime_vec)) %>% 
    mutate(n = n_vec) %>% 
    gather(method, time, -n)
})

names(runtime.result) <- X_types

record <- list(knockoffs = knockoffs.name,
               statistic = statistic.name)

save(runtime.result, record, file = here("data", paste0(experiment, ".RData")))


draw_scale_n_curve(experiment, X_types)





