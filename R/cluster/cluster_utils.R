library(here)

get_expr_index <- function(expr_id, index_data){
    X_chunk_size <- index_data$fig_x_len * index_data$sample_len
    
    X_index <- (expr_id-1) %/% X_chunk_size + 1
    
    fig_x_index <- (expr_id-1) %% X_chunk_size + 1
    fig_x_index <- (fig_x_index-1) %/% index_data$sample_len + 1
    
    y_index <- (expr_id-1) %% index_data$sample_len + 1
    
    return(list(X_index = X_index, fig_x_index = fig_x_index,
                y_index = y_index, expr_id = expr_id))
}

get_data_indeces <- function(job_id, index_data){
    expr_num <- index_data$X_len * index_data$fig_x_len * index_data$sample_len
    exprs_per_job <- ceiling(expr_num / index_data$n_jobs)
    
    job_exprs <- ((job_id-1) * exprs_per_job + 1) : (job_id * exprs_per_job)
    job_exprs <- job_exprs[job_exprs <= expr_num]
    
    job_expr_indeces <- lapply(job_exprs, function(expr_id){
        get_expr_index(expr_id, index_data)
    })
    
    return(job_expr_indeces)
}

