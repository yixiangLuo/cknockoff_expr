library(here)

library(tidyverse)
library(knockoff)
library(cknockoff)

source(here("R", "methods.R"))
source(here("R", "utils.R"))
source(here("R", "plot.R"))


experiment <- "scale_expr"

X_types <- c("IID_Normal", "MCC")
alpha <- 0.05
p_vec <- c(100, 200, 500, 1000, 2000)
pi1_list <- list(fixed_alt = 10 / p_vec, fixed_ratio = rep(0.1, length(p_vec)))
n_rounds <- c(10, 10, 5, 1, 1)

target <- 0.5
target_at_alpha <- 0.2

n_cores <- 7

X_seed <- 2021
noise <- quote(rnorm(n))

knockoffs.name <- quote(create.fixed)
statistic.name <- quote(stat.glmnet_coefdiff_lm)

knockoffs <- get(knockoffs.name)
statistic <- get(statistic.name)

runtime.result <- lapply(X_types, function(X_type){
  print(X_type)
  
  posit_type <- ifelse(X_type != "Homo_Block", "random", "fix")
  
  lapply(pi1_list, function(pi1_vec){
    print(pi1_vec)
    
    runtime_vec <- sapply(1:length(p_vec), function(iter){
      p <- p_vec[iter]
      n <- 3*(p+1)
      
      print(p)
      
      n_round <- n_rounds[iter]
      
      X <- gene_X(X_type, n, p, X_seed)
      
      pi1 <- pi1_vec[iter]
      mu1 <- BH_lm_calib(X, pi1, noise, posit_type, 1, side = "two", nreps = 50,
                         alpha = target_at_alpha, target = target, n_cores = n_cores)
      beta <- genmu(p, pi1, mu1, posit_type, 1)
      
      H0 <- beta == 0
      
      print(paste0(c("alt:", which(!H0)), collapse = " "))
      
      runtime <- sapply(1:n_round, function(iter){
        
        y <- X %*% beta + eval(noise)
        
        time_start <- Sys.time()
        kn.result <- knockoff.filter(X, y, knockoffs = knockoffs,
                                     statistic = statistic, fdr = alpha)
        time_end <- Sys.time()
        
        kn.runtime <- difftime(time_end, time_start, units = "secs")[[1]]
        
        time_start <- Sys.time()
        ckn.result <- cknockoff(X, y,
                                statistic = statistic,
                                alpha = alpha,
                                kappa = 1,
                                n_cores = 1)
        time_end <- Sys.time()
        
        ckn.runtime <- difftime(time_end, time_start, units = "secs")[[1]]
        
        print(paste0(c("kn:", kn.result$selected), collapse = " "))
        print(paste0(c("ckn:", ckn.result$selected), collapse = " "))
        
        runtime_sample <- c(kn.runtime, ckn.runtime)
        names(runtime_sample) <- c("knockoff", "cKnockoff")
        
        return(runtime_sample)
      })
      runtime_avg <- rowMeans(runtime)
      
      return(runtime_avg)
    })
    
    runtime <- as.data.frame(t(runtime_vec)) %>% 
      mutate(p = p_vec) %>% 
      gather(method, time, -p)
  })
})

names(runtime.result) <- X_types

record <- list(pi1_list = pi1_list,
               knockoffs = knockoffs.name,
               statistic = statistic.name)

save(runtime.result, record, file = here("data", paste0(experiment, ".Rdata")))


draw_runtime_curve(experiment, X_types)





