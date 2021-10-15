library(abind)    # abind

library(foreach)
library(doParallel)


# compute fdr and power for methods on a linear problem
get_fdp_power <- function(X, beta, H0, noise = quote(rnorm(n)),
                          alphas, method_list,
                          sample_size = 100, n_cores,
                          expr_name, X_title){
    n <- NROW(X)
    p <- NCOL(X)

    registerDoParallel(n_cores)
    
    org_data <- list()
    FDR_Power <- NULL

    for(alpha_i in 1:length(alphas)){
        alpha <- alphas[alpha_i]
        
        update_progress(expr_name, X_title, alpha, sample_size, action = "start")
        
        results <- foreach(iter = 1:sample_size) %dopar% {
            
            set.seed(iter)
            
            y <- X %*% beta + eval(noise)
            
            sign_beta <- sign(beta)
            
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
            
            update_progress(expr_name, X_title, alpha, sample_size, action = "progress")
            
            return(fdp_power)
        }
        
        update_progress(expr_name, X_title, alpha, sample_size, action = "end")
        
        method_names <- colnames(results[[1]])
        results <- abind(results, along=3)
        
        org_data[[alpha_i]] <- results
        
        interp_results <- function(results, method_names, type, alpha, direction){
            results <- unname(results)
            interpretation <- data.frame(methods = factor(1:NROW(results), ordered = T),
                                         method_name = method_names,
                                         mean = rowMeans(results),
                                         std = sqrt(rowMeans((results - rowMeans(results))^2)),
                                         type = type,
                                         direction = direction,
                                         alpha = alpha)
        }
        
        fdr_power <- rbind(interp_results(results[1, , ], method_names, "FDR", alpha, direction = F),
                           interp_results(results[2, , ], method_names, "Power", alpha, direction = F),
                           interp_results(results[3, , ], method_names, "FDR", alpha, direction = T),
                           interp_results(results[4, , ], method_names, "Power", alpha, direction = T))
        FDR_Power <- rbind(FDR_Power, fdr_power)
    }

    return(list(FDR_Power = FDR_Power, org_data = org_data))
}







