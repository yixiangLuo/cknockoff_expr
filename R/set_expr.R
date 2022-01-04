#!/usr/bin/env Rscript

library(here)

source(here("R", "methods.R"))
source(here("R", "utils.R"))

# source(here("R", "settings", "test.R"))

args <- commandArgs(trailingOnly=TRUE)
if (length(args) == 2) {
  source(here("R", "settings", paste0(args[1], ".R")))
  n_jobs <- as.integer(args[2])
} else {
  stop("accept two paramters: expr_setting and maximal number of jobs.")
}

X_data <- list()
y_data <- list()

for(X_iter in 1:length(X_types)){
  X_type <- X_types[X_iter]
  posit_type <- posit_types[X_iter]
  X_title <- X_titles[X_iter]
  
  X <- gene_X(X_type, n, p, X_seed)
  
  mu1 <- BH_lm_calib(X, pi1, noise = quote(rnorm(n)),
                     posit_type, 1, side = "two", nreps = 200,
                     alpha = target_at_alpha, target = target, n_cores = n_cores)
  beta <- genmu(p, pi1, mu1, posit_type, 1)
  
  H0 <- beta == 0
  
  method_list <- get_method_list(X, knockoffs, statistic)[method_names]
  
  X_data[[X_iter]] <- list(X = X, beta = beta, H0 = H0, method_list = method_list)
  
  y_data[[X_iter]] <- list()
  for(var_i in 1:length(fig_x_var$value)){
    fig_x_value <- fig_x_var$value[var_i]
    alpha <- alphas[[var_i]]
    beta_permute <- beta_permutes[[var_i]]
    noise <- noises[[var_i]]
    
    y_data[[X_iter]][[var_i]] <- list(alpha = alpha, data = list())
    for(iter in 1:sample_size){
      eval(beta_permute)
      y <- X %*% beta + eval(noise)
      
      y_data[[X_iter]][[var_i]]$data[[iter]] <- list(y = y, beta = beta)
    }
  }
}

names(X_data) <- X_types
names(y_data) <- X_types


index_data <- list(X_len = length(X_types), fig_x_len = length(fig_x_var$value),
                   sample_len = sample_size, n_jobs = n_jobs)
expr_num <- length(X_types) * length(fig_x_var$value) * sample_size

unlink(here("data", "temp", experiment), recursive = TRUE) 
dir.create(here("data", "temp", experiment), showWarnings = F)
dir.create(here("data", "temp", experiment, "progress"), showWarnings = F)
save(X_data, y_data, index_data,
     file = here("data", "temp", experiment, "settings.Rdata"))
save(expr_num, file = here("data", "temp", experiment, "expr_num.Rdata"))

unlink(here("log", experiment), recursive = TRUE) 
dir.create(here("log", experiment), showWarnings = F)


## create bash file
# read template
bash_temp <- here("jobs", "expr_template.sh")
con <- file(bash_temp, open="r")
expr_bash <- readLines(con)
close(con)

# specify file content
expr_bash[3] <- paste0("#SBATCH --job-name=", experiment)
expr_bash[4] <- paste0("#SBATCH --output=../log/", experiment, "/%a.out")
expr_bash[5] <- paste0("#SBATCH --array=1-", n_jobs)
expr_bash[14] <- paste0("EXPERIMENT=", experiment)

# write to bash file
expr_bash_file <- here("jobs", paste0(experiment, "_job.sh"))
con <- file(expr_bash_file, open="w")
writeLines(expr_bash, con)
close(con)

