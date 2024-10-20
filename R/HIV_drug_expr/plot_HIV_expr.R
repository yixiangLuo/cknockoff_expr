#### modified from https://github.com/lihualei71/dbhPaper/blob/master/R/plot_dBH_HIV.R

library(here)
library(tidyverse)
library(latex2exp)

source(here("R", "utils.R"))

load(here("data", "HIV", "HIV_discoveries.RData"))

# methods_levels <- c("BH", "dBH", "knockoff", "cKnockoff", "cKnockoff_STAR")
# methods_labels <- parse_name(c("BH", "dBH", "KN", "cKN", "cKN_STAR"))
methods_levels <- c("BH", "knockoff", "cKnockoff", "cKnockoff_STAR")
methods_labels <- parse_name(c("BH", "KN", "cKN", "cKN_STAR"))
names(methods_labels) <- methods_levels

multi_method_color <- c("#984ea3", "dodgerblue3", "#333333", "#33a02c", "red", "orange1")
names(multi_method_color) <- c("BH", "dBH", "knockoff", "BonBH", "cKnockoff", "cKnockoff_STAR")
method_colors <- multi_method_color[methods_levels]

for (al in unique(discoveries$alpha)){
  data <- discoveries %>%
    filter(alpha == al,
           method %in% methods_levels) %>%
    mutate(method = factor(method,
                           levels = methods_levels)) %>%
    mutate(drug = paste0(drug_name, " (", drug_class, ")"))
  drug_levels <- unique(data$drug)
  
  # TSM_num <- data.frame(drug = factor(unique(data$drug), levels = drug_levels),
  #                       m1 = c(rep(34, 7), rep(24, 6), rep(15, 3)))
  TSM_num <- data.frame(drug = factor(unique(data$drug), levels = drug_levels),
                        m1 = c(rep(60, 7), rep(43, 6), rep(26, 3)))
  
  plot <- data %>%
    select(-alpha, -drug_name, -drug_class) %>%
    gather("discoveries", "value", -drug, -method) %>%
    mutate(discoveries = factor(
      discoveries,
      levels = c("nfalse", "ntrue"),
      labels = c("Not in TSM list", "In TSM list")),
      drug = factor(drug,
                    levels = drug_levels)) %>%
    ggplot(aes(x = method, y = value)) +
    geom_bar(stat = "identity", aes(color = method, fill = method, alpha = discoveries)) +
    geom_hline(data = TSM_num, aes(yintercept = m1),
               linetype = "longdash", alpha = 0.6, na.rm = T) +
    facet_wrap(~ drug, nrow = 4) +
    scale_x_discrete(labels = methods_labels) +
    scale_fill_manual(values = method_colors) +
    scale_colour_manual(values = method_colors) +
    guides(fill = "none", color = "none") +
    ylab("Number of discoveries") + 
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(angle = 90),
          strip.text = element_text(size = 10),
          axis.text = element_text(size = 10),
          axis.title = element_text(size = 12.5),
          legend.position = "bottom")
  ggsave(filename = here("figs", paste0("HIV-", al, ".pdf")),
         plot, width = 6.5, height = 7.5)
  
}