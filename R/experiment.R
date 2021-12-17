library(here)

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
load(here("data", "temp", experiment, "settings.Rdata"))

# the indeces of experiments this job should do
job_expr_indeces <- get_data_indeces(job_id, index_data)

# release memory of unused data
for(X_iter in 1:index_data$X_len){
  if(!(X_iter %in% sapply(job_expr_indeces, function(ind){ind$X_index}))){
    X_data[[X_iter]] <- NA
    y_data[[X_iter]] <- NA
  } else{
    for(var_i in 1:index_data$fig_x_len){
      if(!(var_i %in% sapply(job_expr_indeces, function(ind){ind$fig_x_index}))){
        y_data[[X_iter]][[var_i]] <- NA
      } else{
        for(iter in 1:index_data$sample_len){
          if(!(iter %in% sapply(job_expr_indeces, function(ind){ind$y_index}))){
            y_data[[X_iter]][[var_i]]$data[[iter]] <- NA
          }
        }
      }
    }
  }
}

# do the experiments
for(expr_index in job_expr_indeces){
  # record the start of the program
  cat(as.integer(Sys.time()),
      file = here("data", "temp", experiment, "progress", expr_index$expr_id))
  
  X <- X_data[[expr_index$X_index]]$X
  method_list <- X_data[[expr_index$X_index]]$method_list
  
  alpha <- y_data[[expr_index$X_index]][[expr_index$fig_x_index]]$alpha
  y <- y_data[[expr_index$X_index]][[expr_index$fig_x_index]]$data[[expr_index$y_index]]$y
  
  expr_result <- lapply(method_list, function(method){
    method(y, X, alpha)
  })
  
  # save results
  save(expr_result,
       file = here("data", "temp", experiment, paste0(expr_index$expr_id, ".Rdata")))
  
  # delete the "in-progress" record
  unlink(here("data", "temp", experiment, "progress", expr_index$expr_id))
}







