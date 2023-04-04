library(abind)    # abind

library(foreach)
library(doParallel)


# compute fdr and power for methods on a linear problem
get_fdp_power <- function(problem_setting, beta_permutes = NA,
                          noises = quote(rnorm(n)), alphas,
                          fig_x_var,
                          sample_size = 100, n_cores,
                          expr_name, X_title){
    
    mu1 <- problem_setting$mu1
    if(!problem_setting$random_X){
        X <- problem_setting$X
        method_list <- get_method_list(X, problem_setting$knockoffs,
                                       problem_setting$statistic,
                                       problem_setting$method_names)
    }

    registerDoParallel(n_cores)

    org_data <- list()
    FDR_Power <- NULL

    for(var_i in 1:length(fig_x_var$value)){
        fig_x_value <- fig_x_var$value[var_i]
        alpha <- alphas[[var_i]]
        beta_permute <- beta_permutes[[var_i]]
        noise <- noises[[var_i]]

        update_progress(expr_name, X_title, alpha, sample_size, action = "start")

        results <- foreach(iter = 1:sample_size, .options.multicore = list(preschedule = F)) %dopar% {
        # results <- lapply(1:sample_size, function(iter){
            print(iter)

            set.seed(iter)
            
            if(problem_setting$random_X){
                X.res <- gene_X(problem_setting$X_type,
                                problem_setting$n,
                                problem_setting$p,
                                iter)
                X <- X.res$X
                method_list <- get_method_list(X, problem_setting$knockoffs,
                                               problem_setting$statistic,
                                               problem_setting$method_names,
                                               X.res$Xcov.true)
            }
            
            beta <- genmu(problem_setting$p, problem_setting$pi1,
                          mu1, problem_setting$posit_type, 1)
            H0 <- beta == 0
          
            eval(beta_permute)

            # y <- X %*% beta + eval(noise)
            suc_prob <- exp(X %*% beta) / (exp(X %*% beta) + 1)
            y <- sapply(1:n, function(obs_i){
                rbinom(1, 1, suc_prob[obs_i])
            })
            # save(X, y, alpha, file = "debug.RData")

            sign_beta <- sign(beta)

            fdp_power <- NULL
            method_names <- NULL

            for(method_i in 1:length(method_list)){
                # print(names(method_list[method_i]))
                method_res <- method_list[[method_i]](y, X, alpha)
                method_fdp_power <- calc_FDP_power(method_res$selected, H0,
                                                   method_res$sign_predict, sign_beta)
                fdp_power <- cbind(fdp_power, method_fdp_power)
                method_names <- c(method_names, names(method_list[method_i]))
            }
            print(iter)

            colnames(fdp_power) <- method_names

            update_progress(expr_name, X_title, alpha, sample_size, action = "progress")

            return(fdp_power)
        }

        update_progress(expr_name, X_title, alpha, sample_size, action = "end")

        method_names <- colnames(results[[1]])
        results <- abind(results, along=3)

        org_data[[var_i]] <- results

        interp_results <- function(results, method_names, type, alpha, direction){
            results <- unname(results)
            interpretation <- data.frame(methods = factor(1:NROW(results), ordered = T),
                                         method_name = method_names,
                                         mean = rowMeans(results),
                                         std = sqrt(rowMeans((results - rowMeans(results))^2)),
                                         type = type,
                                         direction = direction,
                                         alpha = alpha,
                                         fig_x = fig_x_value)
        }

        fdr_power <- rbind(interp_results(results[1, , ], method_names, "FDR", alpha, direction = F),
                           interp_results(results[2, , ], method_names, "Power", alpha, direction = F),
                           interp_results(results[3, , ], method_names, "FDR", alpha, direction = T),
                           interp_results(results[4, , ], method_names, "Power", alpha, direction = T))
        FDR_Power <- rbind(FDR_Power, fdr_power)
    }

    return(list(FDR_Power = FDR_Power, org_data = org_data))
}







