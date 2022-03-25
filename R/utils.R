# library(mvtnorm)  # rmvnorm
library(here)

library(foreach)
library(doParallel)

# generate design matrix randomly
gene_X <- function(X_type = "IID_Normal", n, p, X_seed = 1){
    cor_radius <- 5
    set.seed(X_seed)
    if(X_type == "IID_Normal"){
        X <- matrix(rnorm(n*p), n) / sqrt(n)
    } else if(X_type == "Coef_AR"){
        rho <- 0.5
        
        cov_mat <- solve(rho^(abs(outer(1:p, 1:p, "-"))))
        
        R <- chol(cov_mat)
        basis <- qr.Q(qr(matrix(rnorm(n*p), n)))
        X <- basis %*% R
    } else if(X_type == "X_AR"){
        rho <- 0.5
        
        cov_mat <- rho^(abs(outer(1:p, 1:p, "-")))
        
        R <- chol(cov_mat)
        basis <- qr.Q(qr(matrix(rnorm(n*p), n)))
        X <- basis %*% R
    } else if(X_type == "Homo_Block"){
        rho <- 0.5
        block_size <- 10

        blockSigma <- matrix(rho, block_size, block_size)
        diag(blockSigma) <- 1

        cov_mat <- as.matrix(diag(p / block_size) %x% blockSigma)
        
        R <- chol(cov_mat)
        basis <- qr.Q(qr(matrix(rnorm(n*p), n)))
        X <- basis %*% R
    } else if(X_type == "MCC"){
        if(n %% (p+1) == 0){
            X <- lapply(1:(n/(p+1)), function(i){
                rbind(diag(rep(1, p)), rep(0, p))
            })
            X <- do.call(rbind, X)
            X <- scale(X, center = T, scale = F)
            X <- scale(X, center = F, scale = sqrt(colSums(X^2)))
        } else{
            cov_mat <- matrix(-1/p, p, p)
            diag(cov_mat) <- 1
            
            R <- chol(cov_mat)
            basis <- qr.Q(qr(matrix(rnorm(n*p), n)))
            X <- basis %*% R
        }
    } else if(X_type == "MCC_Block"){
        block_size <- 5
        
        blockSigma <- matrix(-1/block_size, block_size, block_size)
        diag(blockSigma) <- 1
        
        cov_mat <- as.matrix(diag(p / block_size) %x% blockSigma)
        
        R <- chol(cov_mat)
        basis <- qr.Q(qr(matrix(rnorm(n*p), n)))
        X <- basis %*% R
    } else if(X_type == "Sparse"){
        sparsity <- 0.01
        X <- diag(1, nrow = n, ncol = p)
        lower_tri <- lower.tri(X)
        X[lower_tri] <- replicate(sum(lower_tri), rbinom(1, 1, sparsity))
    }
    # X <- scale(X, center = FALSE, scale = sqrt(colSums(X^2)))

    return(X)
}

corr_noise <- function(n, rho){
    cov_mat <- matrix(0, n, n)
    cov_mat[abs(row(cov_mat) - col(cov_mat)) <= 1] <- -rho * 0.5
    diag(cov_mat) <- 1
    
    R <- chol(cov_mat)
    noise <- t(R) %*% matrix(rnorm(n), nrow = n)
    
    return(noise)
}

makeup_vectors <- function(...){
    vars <- list(...)
    max_length <- max(sapply(vars, length))
    
    env <- parent.frame()
    for(var_i in 1:length(vars)){
        var_name <- names(vars[var_i])
        var_length <- length(vars[[var_i]])
        
        if(var_length < max_length){
            var_list <- list()
            for(i in 1:max_length) var_list <- c(var_list, vars[[var_i]])
            
            env[[var_name]] <- var_list
        } else{
            env[[var_name]] <- as.list(vars[[var_i]])
        }
    }
    
    invisible()
}


calc_FDP_power <- function(rejs, H0, sign_predict = NULL, sign_beta = NULL){
    nrejs <- length(rejs)

    false_discovery <- length(intersect(rejs, which(H0)))
    true_discovery <- nrejs - false_discovery

    FDP <- false_discovery / max(nrejs, 1)
    power <- true_discovery / max(sum(!H0), 1)

    if(!is.null(sign_predict)){
        false_dir <- sum((sign_predict * sign_beta)[rejs] <= 0)
        true_dir <- nrejs - false_dir

        FDP_dir <- false_dir / max(nrejs, 1)
        power_dir <- true_dir / length(sign_beta)
    } else{
        FDP_dir <- NA
        power_dir <- NA
    }

    return(c(FDP, power, FDP_dir, power_dir))
}


## calibrate signal strength, modified from Lihua
BH_lm_calib <- function(X, random_X.data,
                        pi1, noise = quote(rnorm(n)),
                        mu_posit_type, mu_size_type,
                        side,
                        nreps = 200,
                        alpha = 0.05,
                        target = 0.3,
                        n_cores = 7){
    if(!random_X.data$random_X){
        n <- nrow(X)
        p <- ncol(X)
    } else{
        n <- random_X.data$n
        p <- random_X.data$p
    }
    
    if(!random_X.data$random_X){
        X_list <- list(X)
        Sigma_list <- list(solve(t(X) %*% X))
    } else{
        X_list <- lapply(1:random_X.data$sample_num, function(i){
            gene_X(random_X.data$X_type, n, p, i)
        })
        Sigma_list <- lapply(X_list, function(X){
            solve(t(X) %*% X)
        })
    }
    X_sample_num <- length(X_list)
    
    beta_list <- lapply(1:nreps, function(i){
        beta <- genmu(p, pi1, 1, mu_posit_type, mu_size_type)
        if (side == "right"){
            beta <- abs(beta)
        } else if (side == "left"){
            beta <- -abs(beta)
        }
        return(beta)
    })
    eps_list <- lapply(1:nreps, function(i){
        eval(noise)
    })
    
    registerDoParallel(n_cores)
    
    BH_power <- function(mu1){
        power <- unlist(foreach(i = 1:nreps) %dopar% {
            H0 <- beta_list[[i]] == 0            
            beta <- beta_list[[i]] * mu1
            eps <- eps_list[[i]]
            
            X <- X_list[[(i%%X_sample_num)+1]]
            Sigma <- Sigma_list[[(i%%X_sample_num)+1]]
            
            y <- X %*% beta + eps
            
            rejs_BH <- BH_lm(y, X, side = "two", alpha, Sigma = Sigma)$rejs
            power_sample <- calc_FDP_power(rejs_BH, H0)[2]
            
            return(power_sample)
        })
        mean(power) - target
    }
    
    lower <- 0
    upper <- 10
    while (TRUE & upper < 1000){
        tmp <- try(uniroot(BH_power, c(lower, upper))$root)
        if (class(tmp) == "try-error"){
            upper <- upper * 2
        } else {
            return(tmp)
        }
    }
    return(NA)
}

genmu <- function(n, pi1, mu1,
                  posit_type = c("random", "rand_block5", "equi", "head"),
                  mu_type = 1:3){
    m <- ceiling(n * pi1)
    posit_type <- posit_type[1]
    mu_type <- mu_type[1]
    if (posit_type == "random"){
        inds <- sample(n, m, replace = F)
    } else if (posit_type == "rand_block5"){
        block_size <- 5
        if(n %% block_size != 0) stop("#variables not a multiple of block size 5")
        if(m > n/block_size) stop("#non-null is greater than #blocks")
        inds <- (sample(n/block_size, m, replace = F) - 1) * block_size + 1
    } else if (posit_type == "equi"){
        inds <- seq(1, n, floor(1 / pi1))[1:m]
    } else if (posit_type == "head"){
        inds <- 1:m
    }
    mu <- rep(0, n)
    altmu <- switch(mu_type,
                    `1` = rep(1, m),
                    `2` = rnorm(m),
                    `3` = rep(1, m) + 0.15 * (2 * rbinom(m, 1, 0.5) - 1))
    mu[inds] <- mu1 * altmu
    mu
}

try_repeat <- function(expr, default = NULL, n_times = 100){
    success <- F
    for(iter in 1:n_times){
        tryCatch({
            eval(expr, envir = parent.frame())
            success <- T
        }, error = function(msg){}, warning = function(msg){})
        if(success) break
        Sys.sleep(0.01)
    }
    if(!success) eval(default, envir = parent.frame())
    
    return(success)
}

update_count <- function(expr_name, action){
    file_name <- here("data", "temp", paste0("progress-", expr_name, ".RData"))
    
    if(action == "start"){
        iters_done <- 0
        start_time <- Sys.time()
        save(iters_done, start_time, file = file_name)
    } else if(action == "progress"){
        success <- try_repeat(load(file = file_name))
        if(success){
            iters_done <- iters_done + 1
            try_repeat(save(iters_done, start_time, file = file_name))
        } else{
            iters_done <- "(can't access progress record)"
        }
    } else if(action == "end"){
        load(file = file_name)
        iters_done <- Sys.time() - start_time  # abuse of var name
        file.remove(file_name)
    }
    
    return(iters_done)
}

print_progress <- function(expr_name, X_title, alpha, iters_done, sample_size, action){
    file_name <- here("data", "temp", paste0("progress-", expr_name, ".txt"))
    
    try_num <- ifelse(action == "progress", 100, 1)
    try_repeat(progress_record <- readLines(file_name), 
               default = {progress_record <- NULL},
               n_times = try_num)
    
    if(action == "start"){
        progress_record <- c(progress_record,
                             paste0("-------------- ", X_title, " --------------"),
                             paste0("alpha: ", alpha, " -- start: ", Sys.time()),
                             paste0("alpha: ", alpha, " -- 0/", sample_size, " done."))
        writeLines(progress_record, con = file_name)
    } else if(action == "progress"){
        if(length(progress_record) > 1){
            progress_record <- progress_record[1:(length(progress_record)-1)]
            progress_record <- c(progress_record,
                                 paste0("alpha: ", alpha, " -- ", iters_done, "/", sample_size, " done."))
            try_repeat(writeLines(progress_record, con = file_name))
        }
    } else if(action == "end"){
        runtime <- round(as.double(iters_done), digits = 2)
        time_unit <- units(iters_done)
        progress_record <- c(progress_record,
                             paste0("alpha: ", alpha, " -- end: ", Sys.time()),
                             paste("--------- runtime:", runtime, time_unit, "---------"),
                             "")
        writeLines(progress_record, con = file_name)
    }
    
    
}

update_progress <- function(expr_name, X_title, alpha, sample_size, action){
    iters_done <- update_count(expr_name, action)
    print_progress(expr_name, X_title, alpha, iters_done, sample_size, action)
}

# borrowed from knockoff
# Evaluate an expression with the given random seed, then restore the old seed
with_seed = function(seed, expr) {
    seed.old = if (exists('.Random.seed')) .Random.seed else NULL
    set.seed(seed)
    on.exit({
        if (is.null(seed.old)) {
            if (exists('.Random.seed'))
                rm(.Random.seed, envir=.GlobalEnv)
        } else {
            .Random.seed <<- seed.old
        }
    })
    expr
}


# convert a digit to a string
digit_to_char <- function(number){
    return(str_replace(as.character(number), "\\.", "d"))
}

# convert a string to a digit
char_to_digit <- function(char){
    return(str_replace(as.character(char), "d", "\\."))
}

char_to_digit <- function(str){
    gsub(paste0("([0-9]+)d([0-9]+)"),
         paste0("\\1.\\2"),
         str)
}

parse_name <- function(str){
    str <- char_to_digit(str)
    str <- str_replace(str, "_L_", "(")
    str <- str_replace(str, "_R_", ")")
    str <- str_replace(str, "_D_", "-")
    str <- str_replace(str, "_STAR", "*")
}

