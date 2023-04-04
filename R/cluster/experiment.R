library(here)

source(here("R", "cluster", "cluster_utils.R"))

# read cmd args
args <- commandArgs(trailingOnly=TRUE)
if (length(args) == 2) {
  experiment <- args[1]
  job_id <- as.integer(args[2])
} else {
  stop("accept two paramters: expr_name and job id.")
}

# experiment <- "test"
# job_id <- 1

# load packages and data
source(here("R", "settings", paste0(experiment, ".R")))
load(here("data", "temp", experiment, "settings.RData"))

# the indeces of experiments this job should do
job_expr_indeces <- get_data_indeces(job_id, index_data)

# release memory of unused data
for(X_iter in 1:index_data$X_len){
  if(!(X_iter %in% sapply(job_expr_indeces, function(ind){ind$X_index}))){
    X_data[[X_iter]] <- NA
    y_data[[X_iter]] <- NA
  }
}

# do the experiments
for(expr_index in job_expr_indeces){
  # record the start of the program
  cat(as.integer(Sys.time()),
      file = here("data", "temp", experiment, "progress", expr_index$expr_id))
    
  set.seed(expr_index$expr_id)
  
  setting <- X_data[[expr_index$X_index]]
  if(setting$random_X){
      X.res <- gene_X(setting$X_type, setting$n, setting$p, expr_index$expr_id)
      X <- X.res$X
      method_list <- get_method_list(X, setting$knockoffs, setting$statistic,
                                     setting$method_names,
                                     X.res$Xcov.true)
  } else{
      X <- setting$X
      method_list <- setting$method_list
  }
  mu1 <- setting$mu1
  
  beta <- genmu(setting$p, setting$pi1, setting$mu1, setting$posit_type, 1)
  H0 <- beta == 0
  
  alpha <- setting$alphas[[expr_index$fig_x_index]]
  beta_permute <- setting$beta_permutes[[expr_index$fig_x_index]]
  noise <- setting$noises[[expr_index$fig_x_index]]
  
  eval(beta_permute)
  sign_beta <- sign(beta)
  
  y <- X %*% beta + eval(noise)

  fdp_power <- NULL
  method_names <- NULL
  
  for(method_i in 1:length(method_list)){
      method_res <- method_list[[method_i]](y, X, alpha)
      method_fdp_power <- calc_FDP_power(method_res$selected, H0,
                                         method_res$sign_predict, sign_beta)
      fdp_power <- cbind(fdp_power, method_fdp_power)
      method_names <- c(method_names, names(method_list[method_i]))
  }
  
  colnames(fdp_power) <- method_names
  
  # save results
  save(fdp_power,
       file = here("data", "temp", experiment, paste0(expr_index$expr_id, ".RData")))
  
  # delete the "in-progress" record
  unlink(here("data", "temp", experiment, "progress", expr_index$expr_id))
}







