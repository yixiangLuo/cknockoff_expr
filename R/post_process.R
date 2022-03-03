library(here)
library(abind)

source(here("R", "utils.R"))
source(here("R", "cluster_utils.R"))
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

# record computed result in the data structure
expr_num <- index_data$X_len * index_data$fig_x_len * index_data$sample_len
for(expr_id in 1:expr_num){
    expr_index <- get_expr_index(expr_id, index_data)
    
    load(here("data", "temp", experiment, paste0(expr_index$expr_id, ".Rdata")))
    y_data[[expr_index$X_index]][[expr_index$fig_x_index]]$data[[expr_index$y_index]]$fdp_power <- fdp_power
}

X_results <- lapply(1:index_data$X_len, function(X_iter){
    
    org_data <- list()
    FDR_Power <- NULL
    
    for(var_i in 1:index_data$fig_x_len){
        fig_x_value <- fig_x_var$value[var_i]
        alpha <- alphas[[var_i]]
        data <- y_data[[X_iter]][[var_i]]$data
        
        results <- lapply(1:sample_size, function(iter){
            data[[iter]]$fdp_power
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





