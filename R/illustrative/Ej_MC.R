library(here)
library(cknockoff)
library(latex2exp)
library(tidyverse)

source(here("R", "utils.R"))

source(here("R", "settings", "main_expr_10.R"))

iter_i <- 1
X_type <- X_types[iter_i]
posit_type <- posit_types[iter_i]
X_title <- X_titles[iter_i]
random_X <- random_Xs[iter_i]
if(exists("targets")) target <- targets[iter_i]
if(exists("pi1s")) pi1 <- pi1s[iter_i]

print(paste(X_type, pi1, target))

if(!random_X){
  X <- gene_X(X_type, n, p, X_seed)$X
  random_X.data <- list(random_X = random_X)
  mu1 <- signal_calib(calib_method, X, random_X.data, pi1, noise = quote(rnorm(n)),
                      posit_type, 1, side = "two", nreps = 50,
                      alpha = target_at_alpha, target = target, n_cores = n_cores)
} else{
  X.res <- gene_X(X_type, n, p, X_seed)
  X <- NA
  random_X.data <- list(random_X = random_X,
                        n = n, p = p, X_type = X_type, Xcov.true = X.res$Xcov.true,
                        sample_num = 5)
  mu1 <- signal_calib(calib_method, X, random_X.data, pi1, noise = quote(rnorm(n)),
                      posit_type, 1, side = "two", nreps = 50,
                      alpha = target_at_alpha, target = target, n_cores = n_cores)
}

var_i <- 2

alpha <- alphas[[var_i]]
beta_permute <- beta_permutes[[var_i]]
noise <- noises[[var_i]]

seed <- 1000
if(random_X){
  X.res <- gene_X(X_type, n, p, X_seed = seed)
  X <- X.res$X
}
set.seed(seed+1)

beta <- genmu(p, pi1, mu1, posit_type, 1)
H0 <- beta == 0

eval(beta_permute)

y <- X %*% beta + eval(noise)


ckn.data <- cknockoff(X, y, alpha = 0.1, n_cores = 10)
ckn.Ej <- as.data.frame(t(ckn.data))
colnames(ckn.Ej) <- c('j', 'Ej', 'Ej_mc', 'upper', 'lower', 'Ej_n', 'Ej_mc_n')
ckn.Ej$Variable <- c("Non-null", "Null")[(ckn.Ej$j %in% which(H0))+1]
ckn.Ej$Ej <- ckn.Ej$Ej / (ckn.Ej$upper - ckn.Ej$lower)
ckn.Ej$Ej_mc <- ckn.Ej$Ej_mc / (ckn.Ej$upper - ckn.Ej$lower)

ggplot(filter(ckn.Ej, T)) +
  geom_hline(yintercept = 0) +
  geom_line(aes(x = Ej, y = Ej), linetype = "dashed") + 
  geom_point(aes(x = Ej, y = Ej_mc, color = Variable), shape = 1) + 
  theme_bw() +
  theme(aspect.ratio = 1,
        strip.text = element_text(size = 15),
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 10),
        legend.position = "right",
        legend.title=element_text(size=9),
        legend.text=element_text(size=9)) +
  labs(x = TeX("true $\\widetilde{E}_j$"), y = TeX("Monte-Carlo estimated $\\widetilde{E}_j$"))

ggsave(filename = here("figs", paste0("Ej_MC.pdf")),
       width = 5, height = 4)

save(seed, X, y, beta, H0, alpha, ckn.data,
     file = here("data", "temp", paste0("Ej_mc", "-", alpha, "-", seed, ".RData")))
