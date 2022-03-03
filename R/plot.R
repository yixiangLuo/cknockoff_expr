library(here)
library(tidyverse)
library(latex2exp)
# library(cowplot)

source(here("R", "utils.R"))


# fdp and power of the methods
draw_fdp_power_curve <- function(experiment, X_types, sample_size = 1,
                                 method_names, method_colors, method_shapes,
                                 error_bar = F, direction = F){
    load(here("data", paste0(experiment, ".Rdata")))
    
    fdr_power <- lapply(X_types, function(X_type){
        results[[X_type]]$FDR_Power %>% mutate(design_mat = str_replace(X_type, "_", "-"))
    })
    fdr_power <- do.call(rbind, fdr_power)
    X_types <- str_replace(X_types, "_", "-")
    rm(results)
    
    dir <- direction
    fdr_power <- fdr_power %>% 
        filter(direction == dir,
               design_mat %in% X_types,
               method_name %in% method_names) %>%
        arrange(methods)

    methods_level <- unique(fdr_power$methods)
    # alphas <- unique(fdr_power$alpha)

    method_names <- parse_name(method_names)

    ref_prototype <- data.frame(fig_x = c(fig_x_var$value, fig_x_var$value),
                                threshold = c(unlist(alphas), rep(NA, length(fig_x_var$value))),
                                type = rep(c("FDR", "Power"), each = length(fig_x_var$value)))
    reference <- lapply(X_types, function(X_type){
        ref_prototype %>% mutate(design_mat = X_type)
    })
    reference <- do.call(rbind, reference)
    
    if(error_bar){
        add_error_bar <- quote(
            geom_errorbar(aes(x = fig_x, y = mean,
                              ymin = mean-2*std/sqrt(sample_size),
                              ymax = mean+2*std/sqrt(sample_size),
                              color = methods), width=0.05,
                          position = position_dodge(width=0.01))
        )
    } else{
        add_error_bar <- NULL
    }
    
    if(fig_x_var$value[1] <= fig_x_var$value[length(fig_x_var$value)]){
        set_x_axis <- scale_x_continuous
    } else{
        set_x_axis <- scale_x_reverse
    }

    plot <- ggplot(fdr_power) +
        geom_line(aes(x = fig_x, y = mean, color = methods)) +
        geom_point(aes(x = fig_x, y = mean, color = methods, shape = methods), size = 2) +
        eval(add_error_bar) +
        geom_line(data = reference, aes(x = fig_x, y = threshold),
                  linetype = "longdash", alpha = 0.6, na.rm = T) +
        facet_grid(vars(factor(design_mat, levels = X_types)),
                   vars(factor(type, levels = c("FDR", "Power"))), scales="free") +
        # facet_wrap(vars(factor(design_mat, levels = X_types), factor(type, levels = c("FDR", "Power"))),
        #            ncol = 2, scales="free_y") +
        set_x_axis(breaks = fig_x_var$value, labels = fig_x_var$value) +
        scale_color_manual(values = method_colors, labels = method_names, breaks = methods_level) +
        scale_shape_manual(values = method_shapes, labels = method_names, breaks = methods_level) +
        theme_bw() +
        theme(aspect.ratio = 1,
              panel.grid = element_blank(),
              strip.text = element_text(size = 15),
              axis.title = element_text(size = 13),
              axis.text = element_text(size = 10),
              legend.position = "right",
              legend.title=element_text(size=9),
              legend.text=element_text(size=9)) +
        labs(x = fig_x_var$name, y = "Estimated FDR/Power")

    ggsave(filename = here("figs", paste0("simu-", experiment, ".pdf")),
           plot, width = 7, height = 2*(length(X_types)+1))
    
    return(plot)
}

draw_fdp_power_dist <- function(experiment, X_types, method_names, method_colors,
                                type = "TPP", figure = "ECDF"){
    load(here("data", paste0(experiment, ".Rdata")))
    
    method_names_org <- unique(results[[1]]$FDR_Power$method_name)
    method_level_org <- factor(seq_along(method_names_org), ordered = T)
    names(method_level_org) <- method_names_org
    data_index <- ifelse(type == "FDP", 1, 2)
    cdf_x <- seq(from = 0, to = 1, length.out = 100)
    
    plot_data <- lapply(X_types, function(X_type){
        org_data_list <- results[[X_type]]$org_data
        org_df <- lapply(1:length(org_data_list), function(fig_x_ind){
            data_fig_x <- org_data_list[[fig_x_ind]][data_index, , ]
            rownames(data_fig_x) <- method_names_org
            
            if(figure == "hist"){
                data_fig_x <- as.data.frame(t(data_fig_x)) %>%
                    pivot_longer(everything(), names_to = "method_name", values_to = "value") %>%
                    mutate(fig_x = fig_x_var$value[fig_x_ind],
                           methods = unname(method_level_org[method_name]))
            } else{
                data_fig_x <- sapply(as.data.frame(t(data_fig_x)), function(col){
                    ECDF <- ecdf(col)(cdf_x)
                })
                data_fig_x <- as.data.frame(data_fig_x) %>% 
                    mutate(ecdf_x = cdf_x) %>%
                    pivot_longer(!ecdf_x, names_to = "method_name", values_to = "value") %>%
                    mutate(fig_x = fig_x_var$value[fig_x_ind],
                           methods = unname(method_level_org[method_name]))
            }
            
            return(data_fig_x)
        })
        do.call(rbind, org_df) %>% mutate(design_mat = str_replace(X_type, "_", "-"))
    })
    plot_data <- do.call(rbind, plot_data)
    
    X_types <- str_replace(X_types, "_", "-")
    rm(results)
    
    plot_data <- plot_data %>% 
        filter(design_mat %in% X_types,
               method_name %in% method_names) %>%
        arrange(methods)
    
    methods_level <- unique(plot_data$methods)
    
    method_names <- parse_name(method_names)
    
    if(figure == "hist"){
        plot <- ggplot(plot_data) +
            geom_histogram(aes(x = value, fill = methods, color = methods),
                           alpha = 0.5, position = "identity") +
            facet_grid(vars(factor(design_mat, levels = X_types)),
                       vars(factor(fig_x, levels = fig_x_var$value)), scales="free") +
            scale_fill_manual(values = method_colors, labels = method_names, breaks = methods_level) +
            scale_color_manual(values = method_colors, labels = method_names, breaks = methods_level) +
            theme_bw() +
            theme(aspect.ratio = 1,
                  # panel.grid = element_blank(),
                  strip.text = element_text(size = 13),
                  axis.title = element_text(size = 11),
                  axis.text = element_text(size = 8),
                  legend.position = "right",
                  legend.title=element_text(size=9),
                  legend.text=element_text(size=9)) +
            labs(x = type)
    } else{
        plot <- ggplot(plot_data) +
            geom_line(aes(x = ecdf_x, y = value, color = methods)) +
            facet_grid(vars(factor(design_mat, levels = X_types)),
                       vars(factor(fig_x, levels = fig_x_var$value)), scales="free") +
            scale_color_manual(values = method_colors, labels = method_names, breaks = methods_level) +
            theme_bw() +
            theme(aspect.ratio = 1,
                  # panel.grid = element_blank(),
                  strip.text = element_text(size = 13),
                  axis.title = element_text(size = 11),
                  axis.text = element_text(size = 8),
                  legend.position = "right",
                  legend.title=element_text(size=9),
                  legend.text=element_text(size=9)) +
            labs(x = type, y = "ECDF")
    }
    

    ggsave(filename = here("figs", paste0(experiment, "-", type, "-", figure ,".pdf")),
           plot, width = 11, height = 2*(length(X_types)+1))
    
    return(plot)
}


draw_scale_m_curve <- function(experiment, X_types){
    load(here("data", paste0(experiment, ".Rdata")))
    alt_types <- c("fixed_alt", "fixed_ratio")
    runtime <- lapply(X_types, function(X_type){
        fixed_alt <- runtime.result[[X_type]]$fixed_alt %>% mutate(alt = alt_types[1])
        fixed_ratio <- runtime.result[[X_type]]$fixed_ratio %>% mutate(alt = alt_types[2])
        rbind(fixed_alt, fixed_ratio) %>% mutate(design_mat = str_replace(X_type, "_", "-"))
    })
    runtime <- do.call(rbind, runtime)
    X_types <- str_replace(X_types, "_", "-")
    rm(runtime.result, record)
    
    p_vec <- unique(runtime$p)
    method_names <- unique(runtime$method)
    method_colors <- unname(multi_method_color[method_names])
    method_shapes <- unname(multi_method_shape[method_names])
    method_names <- parse_name(method_names)
    runtime$method <- parse_name(runtime$method)
    
    runtime$design_mat <- factor(runtime$design_mat)
    runtime$alt <- factor(runtime$alt)
    levels(runtime$alt) <- c(fixed_alt = TeX("$#alts = 10$"), fixed_ratio = TeX("$\\pi_1 = 0.1$"))
    
    # temp <- sort(unique(round(runtime$time)), decreasing = T)
    # time_ticks <- temp[1]
    # for(ind in 1:length(temp)){
    #     if(temp[ind] < min(time_ticks) * 0.5)
    #         time_ticks <- c(time_ticks, temp[ind])
    # }
    time_ticks <- seq(0, ceiling(log(max(runtime$time), base = 4)), by = 1)
    time_ticks <- 4^time_ticks
    
    plot <- ggplot(runtime) +
        geom_line(aes(x = p, y = time, color = method)) +
        geom_point(aes(x = p, y = time, color = method, shape = method), size = 2) +
        facet_grid(vars(design_mat), vars(alt),
                   labeller=label_parsed) +
        # scale_x_continuous(breaks = p_vec, labels = p_vec) +
        scale_x_log10(breaks = p_vec, labels = p_vec) +
        scale_y_log10(breaks = time_ticks, labels = time_ticks) +
        scale_color_manual(values = method_colors, labels = method_names, breaks = method_names) +
        scale_shape_manual(values = method_shapes, labels = method_names, breaks = method_names) +
        theme_bw() +
        theme(aspect.ratio = 1,
              # panel.grid = element_blank(),
              strip.text = element_text(size = 13),
              axis.title = element_text(size = 11),
              axis.text = element_text(size = 8),
              legend.position = "right",
              legend.title=element_text(size=9),
              legend.text=element_text(size=9)) +
        labs(x = "Number of Hypotheses: m", y = "Runtime (s)")
    
    ggsave(filename = here("figs", paste0("simu-", experiment, ".pdf")),
           plot, width = 6, height = 5)
    
    return(plot)
}

draw_scale_n_curve <- function(experiment, X_types){
    load(here("data", paste0(experiment, ".Rdata")))
    runtime <- lapply(X_types, function(X_type){
        runtime.result[[X_type]] %>% mutate(design_mat = str_replace(X_type, "_", "-"))
    })
    runtime <- do.call(rbind, runtime)
    X_types <- str_replace(X_types, "_", "-")
    rm(runtime.result, record)
    
    n_vec <- unique(runtime$n)
    method_names <- unique(runtime$method)
    method_colors <- unname(multi_method_color[method_names])
    method_shapes <- unname(multi_method_shape[method_names])
    method_names <- parse_name(method_names)
    runtime$method <- parse_name(runtime$method)
    
    runtime$design_mat <- factor(runtime$design_mat)
    
    # temp <- sort(unique(round(runtime$time)), decreasing = T)
    # time_ticks <- temp[1]
    # for(ind in 1:length(temp)){
    #     if(temp[ind] < min(time_ticks) * 0.5)
    #         time_ticks <- c(time_ticks, temp[ind])
    # }
    time_ticks <- seq(0, ceiling(log(max(runtime$time), base = 4)), by = 1)
    time_ticks <- 4^time_ticks
    
    plot <- ggplot(runtime) +
        geom_line(aes(x = n, y = time, color = method)) +
        geom_point(aes(x = n, y = time, color = method, shape = method), size = 2) +
        facet_grid(NULL, vars(design_mat)) +
        # scale_x_continuous(breaks = n_vec, labels = n_vec) +
        scale_x_log10(breaks = n_vec, labels = n_vec) +
        scale_y_log10(breaks = time_ticks, labels = time_ticks) +
        scale_color_manual(values = method_colors, labels = method_names, breaks = method_names) +
        scale_shape_manual(values = method_shapes, labels = method_names, breaks = method_names) +
        theme_bw() +
        theme(aspect.ratio = 1,
              # panel.grid = element_blank(),
              strip.text = element_text(size = 13),
              axis.title = element_text(size = 11),
              axis.text = element_text(size = 8),
              legend.position = "right",
              legend.title=element_text(size=9),
              legend.text=element_text(size=9)) +
        labs(x = "Number of data points: n", y = "Runtime (s)")
    
    ggsave(filename = here("figs", paste0("simu-", experiment, ".pdf")),
           plot, width = 6, height = 5)
    
    return(plot)
}

draw_kn_conservative <- function(experiment){
    load(here("data", paste0(experiment, ".Rdata")))
    
    quant_names <- rownames(results)[1:(NROW(results)-1)]
    quant_colors <- c("#a6bddb", "#fc8d59", "#2b8cbe", "#d95f0e", "#045a8d", "#333333",
                      "#2b8cbe", "#045a8d", "#333333")
    quant_labels <- unname(TeX(c("$M_0 \\alpha$", "$M_{\\tau} \\alpha$", "$b \\, (H_{0})$", "$M_{\\tau_{0}} \\alpha$", "$b^0 \\, (H_{0})$", "FDR",
                                 "$b \\, (H_{1})$", "$b^0 \\, (H_{1})$", "TDR")))
    solid <- c(0, 0, 1, 0, 0, 0)
    
    quant_to_show <- list(c("M_0", "M_tau0", "b0_F", "FDR", "alpha"),
                          c("M_0", "M_tau", "b_F", "M_tau0", "b0_F", "FDR", "alpha"),
                          c("b_T", "b0_T", "TDR", "alpha"))
    figure_appendix <- c("init", "full", "power")
    
    for(fig_i in 1:length(figure_appendix)){
        kn.conservative <- data.frame(t(results)) %>%
            # mutate(M_0 = pmin(M_0, 1), M_tau = pmin(M_tau, M_0)) %>%
            select(all_of(quant_to_show[[fig_i]])) %>%
            gather(quantity, value, -alpha) %>%
            mutate(quantity = factor(quantity, levels = quant_names, ordered = T))
        
        
        plot <- ggplot(kn.conservative) +
            geom_ribbon(aes(x = alpha, ymax = value, fill = quantity), ymin = 0, alpha = 1) +
            geom_line(aes(x = alpha, y = value, color = quantity, alpha = quantity)) +
            scale_color_manual(values = quant_colors, labels = quant_labels, breaks = quant_names) +
            scale_fill_manual(values = quant_colors, labels = quant_labels, breaks = quant_names) +
            scale_alpha_manual(values = solid, labels = quant_labels, breaks = quant_names) +
            theme_light() +
            theme(aspect.ratio = 1,
                  panel.grid = element_blank(),
                  axis.title = element_text(size = 8),
                  axis.text = element_text(size = 6),
                  legend.position = "right",
                  legend.title = element_text(size=8),
                  legend.text = element_text(size=7),
                  legend.key.size = unit(0.5, 'cm')) +
            labs(x = TeX("$nominal \\, FDR \\, level: \\, \\alpha$"), y = TeX("$value / \\alpha$"))
        
        ggsave(filename = here("figs", paste0("kn_conservative-", figure_appendix[fig_i], ".pdf")),
               plot, width = 4, height = 3)
    }
    
    invisible()
}




