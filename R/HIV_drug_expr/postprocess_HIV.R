#### modified from https://github.com/lihualei71/dbhPaper/blob/master/R/postprocess_dBH_HIV.R

library(here)

load(here("data", "HIV", "HIV_res.RData"))
load(here("data", "HIV", "HIV_data.RData"))

get_position <- function(x)
  sapply(regmatches(x, regexpr('[[:digit:]]+', x)), as.numeric)

## Getting drug names and classes
drug_names <- sapply(data, function(dat){
  names(dat$Y)
})
drug_classes <- sapply(1:length(drug_names), function(i){
  rep(names(drug_names)[i], length(drug_names[[i]]))
})
drug_names <- do.call(c, drug_names)
drug_classes <- do.call(c, drug_classes)

## Getting the positions of the gene list and the TSM gene list for each drug
signal_genes <- list()
for (drug_class in c("PI", "NRTI", "NNRTI")){
  base_url <- 'http://hivdb.stanford.edu/pages/published_analysis/genophenoPNAS2006'
  tsm_url <- paste(base_url, 'MUTATIONLISTS', 'NP_TSM', drug_class, sep = '/')
  tsm_df <- read.delim(tsm_url, header = FALSE, stringsAsFactors = FALSE)
  # signal_genes[[drug_class]] <- rep(list(tsm_df[, 1]),
  #                                   length(data[[drug_class]]$Y))
  signals <- unlist(sapply(1:NROW(tsm_df), function(row){
    paste(tsm_df[row, 1], unlist(strsplit(tsm_df[row, 2], " ")), sep = ".")
  }))
  signal_genes[[drug_class]] <- rep(list(paste0("P", signals)),
                                    length(data[[drug_class]]$Y))
}
signal_genes <- do.call(c, signal_genes)

## Merge results for three types of drug, resulting in a list of length 16
res <- do.call(c, res)

## Replace the rejection indices into the gene position
discoveries <- data.frame()
for (i in 1:length(res)){
  signals <- signal_genes[[i]]
  for (j in 1:length(alphas)){
    num_discoveries <- lapply(1:n_sample, function(k){
      sapply(res[[i]][[k]][[j]]$rejs, function(rejs){
        # rejs <- unique(get_position(rejs))
        # true_discoveries <- length(intersect(rejs, signals))
        # false_discoveries <- length(setdiff(rejs, signals))
        true_discoveries <- sum(rejs %in% signals)
        false_discoveries <- sum(!(rejs %in% signals))
        c(true_discoveries, false_discoveries)
      })
    })
    num_discoveries <- Reduce('+', num_discoveries) / n_sample
    df <- data.frame(
      alpha = res[[i]][[1]][[j]]$alpha,
      drug_class = drug_classes[i],
      drug_name = drug_names[i],
      method = colnames(num_discoveries),
      ntrue = num_discoveries[1, ],
      nfalse = num_discoveries[2, ]
    )
    discoveries <- rbind(discoveries, df)
  }
}

save(res, discoveries, file = here("data", "HIV", "HIV_discoveries.RData"))



