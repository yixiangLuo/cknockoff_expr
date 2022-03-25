#### modified from https://github.com/lihualei71/dbhPaper/blob/master/R/dBH_lm_HIV_expr.R

library(here)
library(dbh)
library(knockoff)
library(cknockoff)

library(glmnet)   # lasso
library(KernSmooth)

source(here("R", "methods.R"))
source(here("R", "utils.R"))

if (!file.exists(here("data", "HIV", "HIV_data.Rdata"))){
  source(here("R", "HIV_drug_expr", "HIV_preprocess.R"))
}
load(here("data", "HIV", "HIV_data.Rdata"))

HIV_expr <- function(X, y, alphas,
                     knockoffs = ckn.create.fixed,
                     statistic = stat.glmnet_coefdiff_tiebreak){
  ## Log-transform the drug resistance measurements.
  y <- log(y)
  
  ## Remove patients with missing measurements.
  missing <- is.na(y)
  y <- y[!missing]
  X <- X[!missing,]
  
  ## Remove predictors that appear less than 3 times.
  X <- X[,colSums(X) >= 3]
  
  ## Remove duplicate predictors.
  X <- X[,colSums(abs(cor(X)-1) < 1e-4) == 1]
  
  ## Get names
  genes <- colnames(X)
  
  ## Get stats
  nalphas <- length(alphas)
  n <- nrow(X)
  p <- ncol(X)
  
  ## process data
  Sigma <- solve(t(X) %*% X)
  X.pack <- process_X(X, knockoffs = knockoffs, intercept = T)
  knockoffs_gene <- function(X){return(X.pack$X_kn.org)}
  
  methods <- c("BH", "dBH", "knockoff", "cKnockoff", "cKnockoff_STAR")
  
  res <- list()
  for (k in 1:nalphas){
    alpha <- alphas[k]
    res[[k]] <- list(alpha = alpha, rejs = list())
    obj <- list()
    
    # BH
    t_result <- lm_to_t(y, X, Sigma)
    pvals <- pvals_t(t_result$tvals, t_result$df, side = "two")
    BH.result <-  BH_weighted(pvals, alpha)
    obj <- c(obj, list(genes[BH.result$rejs]))
    
    # dBH
    dBH.result <- dBH_lm(y, X, intercept = FALSE, side = "two",
                         alpha = alpha, gamma = 0.9, niter = 1,
                         avals_type = "BH", qcap = 2)
    obj <- c(obj, list(genes[dBH.result$rejs]))
    
    # browser()
    # D <- diag(t(X.pack$X) %*% X.pack$X_kn - t(X.pack$X) %*% X.pack$X)
    # View(sort(D))
    
    # knockoff
    y.data <- cknockoff:::transform_y(X.pack, y, intercept = T)
    kn.result <- knockoff.filter(X.pack$X.org, y.data$y.org, knockoffs = knockoffs_gene,
                                 statistic = statistic, fdr = alpha)
    obj <- c(obj, list(genes[kn.result$selected]))
    
    # cKnockoff
    ckn.result <- cknockoff(prelim_result = kn.result,
                            n_cores = 14, X.pack = X.pack)
    obj <- c(obj, list(genes[ckn.result$selected]))
    
    # cKnockoff*
    ckn_star.result <- cknockoff(prelim_result = ckn.result, Rstar_refine = T,
                                 n_cores = 14, X.pack = X.pack)
    obj <- c(obj, list(genes[ckn_star.result$selected]))
    
    names(obj) <- methods
    res[[k]]$rejs <- obj
  }
  
  return(res)
}

set.seed(2021)
n_sample <- 20
alphas <- c(0.05, 0.2)
res <- list()

for (drug_class in names(data)){
  Y <- data[[drug_class]]$Y
  X <- data[[drug_class]]$X
  res[[drug_class]] <- lapply(1:length(Y), function(j){
    print(j)
    lapply(1:n_sample, function(iter){
      HIV_expr(X, Y[[j]], alphas)
    })
  })
  save(res, file = here("data", "HIV", paste0("HIV_res-", drug_class, ".RData")))
}

save(res, n_sample, alphas, file = here("data", "HIV", "HIV_res.RData"))




