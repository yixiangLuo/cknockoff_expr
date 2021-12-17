library(here)
library(abind)

source(here("R", "utils.R"))
source(here("R", "plot.R"))

# read cmd args
args <- commandArgs(trailingOnly=TRUE)
if (length(args) == 1) {
    experiment <- args[1]
} else {
    stop("accept one paramter: expr_name.")
}

# load packages and data
source(here("R", "settings", paste0(experiment, ".R")))
load(here("data", "temp", experiment, "settings.Rdata"))

expr_num <- index_data$X_len * index_data$fig_x_len * index_data$sample_len
for(expr_id in 1:expr_num){
    expr_index <- get_expr_index(expr_id, index_data)
    
    load(here("data", "temp", experiment, paste0(expr_index$expr_id, ".Rdata")))
    y_data[[expr_index$X_index]][[expr_index$fig_x_index]]$data[[expr_index$y_index]]$result <- expr_result
}

X_results <- lapply(1:index_data$X_len, function(X_iter){
    H0 <- X_data[[X_iter]]$H0
    
    org_data <- list()
    FDR_Power <- NULL
    
    for(var_i in 1:index_data$fig_x_len){
        fig_x_value <- fig_x_var$value[var_i]
        alpha <- y_data[[X_iter]][[var_i]]$alpha
        data <- y_data[[X_iter]][[var_i]]$data
        
        results <- lapply(1:sample_size, function(iter){
            sign_beta <- sign(data[[iter]]$beta)
            
            fdp_power <- NULL
            method_names <- NULL
            
            for(method_res in data[[iter]]$result){
                method_fdp_power <- calc_FDP_power(method_res$selected, H0,
                                                   method_res$sign_predict, sign_beta)
                fdp_power <- cbind(fdp_power, method_fdp_power)
            }
            
            colnames(fdp_power) <- method_names
            
            return(fdp_power)
        })
        
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
    
    X_result <- list(FDR_Power = FDR_Power, org_data = org_data)
    
    return(X_result)
})

results <- X_results
names(results) <- X_types

save(results, alphas, fig_x_var, file = here("data", paste0(experiment, ".Rdata")))

draw_fdp_power_curve(experiment, X_types, sample_size,
                     method_names, method_colors, method_shapes,
                     error_bar = F, direction = F)





