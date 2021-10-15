library(here)
library(tidyverse)
library(latex2exp)
# library(cowplot)

source(here("R", "utils.R"))

# fdp and power of the methods
draw_fdp_power <- function(file_name, alpha = 0.2, sample_size = 1,
                           y_lim = NULL, excepts = NULL,
                           ref_colors = c("dodgerblue3", "orange1"),
                           Refadj = F){

    if(Refadj){
        location <- str_locate(file_name, "fdp_power")
        file_name <- insert_to_str(file_name, location[2], "_Refadj")
    }

    fdp_power <- read.csv(paste(file_name, ".csv", sep = "")) %>%
        mutate(methods = factor(methods))

    if(!is.null(excepts)){
        # excepts <- sapply(excepts, function(str){paste0(str, "-")})
        # excepts <- paste0(excepts, collapse = "|")
        fdp_power <- fdp_power %>% filter(!grepl(excepts, method_name))
    }

    way_weight_names <- filter(fdp_power, grepl("-", method_name))$method_name %>% unique()
    way_num <- sub("-.*", "", way_weight_names) %>% unique() %>% length()
    weight_scheme_num <- length(way_weight_names) / max(way_num, 1)

    methods_level <- filter(fdp_power, type == unique(fdp_power$type)[1])$methods
    methods_name <- filter(fdp_power, type == unique(fdp_power$type)[1])$method_name

    color_types <- c("indianred", "lightblue", "lightgoldenrod", "mediumpurple", "palegreen")
    colors <- ref_colors
    if(weight_scheme_num > 0){
        for(scheme_i in 1:weight_scheme_num){
            for(way_i in 1:way_num){
                colors <- c(colors, paste0(color_types[scheme_i], way_i))
            }
        }
    }

    max_ref_power <- fdp_power %>% filter(type == "Power", as.numeric(methods) < 3)
    max_ref_power <- max(max_ref_power$mean)
    reference <- data.frame(type = unique(fdp_power$type), value = c(alpha, max_ref_power))

    if(is.null(y_lim)){
        y_lim <- fdp_power %>% filter(type == "Power")
        y_lim <- max(y_lim$mean)
        y_lim <- max(0.7, ceiling(y_lim * 10) / 10)
    }

    if(length(methods_level) > 7){
        facet_row <- vars(type)
        facet_col <- NULL
    } else{
        facet_row <- NULL
        facet_col <- vars(type)
    }

    ggplot(fdp_power) +
        geom_bar(aes(x=methods, y=mean, fill = methods), stat="identity") +
        geom_errorbar(aes(x=methods, y=mean,
                          ymin=mean-2*std/sqrt(sample_size),
                          ymax=mean+2*std/sqrt(sample_size)), width=0.2) +
        geom_hline(data = reference, aes(yintercept = value),
                   linetype = "longdash", alpha = 0.6) +
        facet_grid(facet_row, facet_col) +
        scale_x_discrete(breaks = methods_level, labels = methods_name) +
        scale_y_continuous(breaks = seq(0, 1, length.out = 6)) +
        scale_fill_manual(values = colors) +
        theme_bw() +
        ylim(0, y_lim) +
        theme(panel.grid = element_blank(),
              axis.text.x = element_text(angle = 60),
              strip.text = element_text(size = 17.5),
              axis.text = element_text(size = 12.5),
              axis.title = element_text(size = 17.5),
              legend.position = "none",
              axis.title.x=element_blank(),
              axis.title.y=element_blank())

    # ggsave(file=paste(file_name, ".eps", sep = ""), width = NA, height = 10, units = "cm")
}

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
    fdr_power <- fdr_power %>% filter(direction == dir,
                                      design_mat %in% X_types,
                                      method_name %in% method_names)

    methods_level <- unique(fdr_power$methods)
    alphas <- unique(fdr_power$alpha)

    method_names <- parse_name(method_names)

    ref_prototype <- data.frame(alpha = c(alphas, alphas),
                                threshold = c(alphas, rep(NA, length(alphas))),
                                type = rep(c("FDR", "Power"), each = length(alphas)))
    reference <- lapply(X_types, function(X_type){
        ref_prototype %>% mutate(design_mat = X_type)
    })
    reference <- do.call(rbind, reference)
    
    if(error_bar){
        add_error_bar <- quote(
            geom_errorbar(aes(x = alpha, y = mean,
                              ymin = mean-2*std/sqrt(sample_size),
                              ymax = mean+2*std/sqrt(sample_size),
                              color = methods), width=0.05,
                          position = position_dodge(width=0.01))
        )
    } else{
        add_error_bar <- NULL
    }

    plot <- ggplot(fdr_power) +
        geom_line(aes(x = alpha, y = mean, color = methods)) +
        geom_point(aes(x = alpha, y = mean, color = methods, shape = methods), size = 2) +
        eval(add_error_bar) +
        geom_line(data = reference, aes(x = alpha, y = threshold),
                  linetype = "longdash", alpha = 0.6, na.rm = T) +
        facet_grid(vars(factor(design_mat, levels = X_types)),
                   vars(factor(type, levels = c("FDR", "Power")))) +
        scale_x_continuous(breaks = alphas, labels = alphas) +
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
        labs(x = "Nominated FDR level", y = "Estimated FDR/Power")

    ggsave(filename = here("figs", paste0("simu-", experiment, ".pdf")),
           plot, width = 7, height = 10)
    
    return(plot)
}

draw_runtime_curve <- function(experiment, X_types){
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




