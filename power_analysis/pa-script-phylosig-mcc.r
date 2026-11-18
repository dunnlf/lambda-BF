### power-analysis script

suppressWarnings({library(truncnorm)
library(optparse)
library(phytools)
library(numDeriv)
library(phangorn)
library(geiger)
library(expm)
library(ape)

})

source("power_analysis/phylosig/phylosig.R")

get_synthetic_data <- function(lam, tree, a=0, sig2=1) {
    C <- get_pagel_cov(lam, tree)

    S <- sqrtm(C)

    x <- rnorm(n=length(tree$tip.label))

    as.numeric(a + sqrt(sig2) * S %*% x)
}

run_power_analysis <- function(trees, lam, tau=NULL, decision_boundary=0.5, N_draws=100, N_samples=25, burnin=0, downsample=1) {
    # decision boundary on log10 bayes factor
    # N_draws: number of data sets to simulate
    # tau: an optional single tree on which data will be simulated. If NULL, a separate tree will be sampled from the input trees at each simulation

    # downsample trees
    if (burnin > 0 || downsample > 1) {
        trees <- trees[seq(burnin + 1, length(trees), downsample)]
    }
    cat(paste("Using", length(trees), "trees after burnin and downsampling"))

    decisions <- c()

    lam_samples <- list()
    p_samples <- list()
    tree_indices <- c()

    for (i in 1:N_draws) {
        # sample a tree
        tree_idx <- sample(seq_along(trees), 1)[[1]]
        tree_indices <- c(tree_indices, tree_idx)

        tree <- trees[[tree_idx]]
        trees_ <- trees[-tree_idx]
        mcc_i <- phangorn::maxCladeCred(trees_)

        x <- setNames(get_synthetic_data(lam, tree), tree$tip.label)
	
        res_lam <- c()
        res_p <- c()

        res <- phylosig(mcc_i, x, method='lambda', test=TRUE)
        res_lam <- c(res_lam, res$lambda)
        res_p <- c(res_p, res$P)
        
        lam_samples[[i]] <- res_lam
        p_samples[[i]] <- res_p

        # make decision based on decision boundary
        decision <- ifelse(mean(res_p) < decision_boundary, 1, 0)
        decisions <- c(decisions, decision)
    }

    lam_samples <- do.call(rbind, lam_samples)
    p_samples <- do.call(rbind, p_samples)

    return(list(samples=lam_samples, p_samples=p_samples, trees=tree_indices, decisions=decisions, power=mean(decisions)))
}


opts <- parse_args(OptionParser(option_list = list(
  make_option("--trees", type = "character"),
  make_option("--lam", type = "double", default = 1.0),
  make_option("--decision_boundary", type = "double", default = 0.05),
  make_option("--N_draws", type = "integer", default = 100),
  make_option("--N_samples", type = "integer", default = 25),
  make_option("--burnin", type = "integer", default = 0),
  make_option("--downsample", type = "integer", default = 1),
  make_option("--output", type = "character", default = NULL),
  make_option("--verbose", type = "logical", default = FALSE)
)))

if (opts$verbose) {
    cat('Running power analysis with \ntree file:',opts$trees,
    '\nlam:', opts$lam,
    '\ndecision boundary:', opts$decision_boundary,
    '\nN draws:', opts$N_draws,
    '\nN samples:', opts$N_samples,
    '\nburnin:', opts$burnin,
    '\ndownsample:', opts$downsample, "\n")
}

if (is.null(opts$trees)) {  # for debugging only
    warning("Debug mode - if you didn't intend this, provide a tree file in nexus format using --trees")
    trees <- lapply(1:10, function(i) rcoal(100))
} else {
    trees <- read.nexus(opts$trees)
}

res <- run_power_analysis(trees, lam=opts$lam, decision_boundary=opts$decision_boundary, N_draws=opts$N_draws,
        N_samples=opts$N_samples, burnin=opts$burnin, downsample=opts$downsample)

tree_file_name <- gsub(".*/", "", opts$trees)

opts$output <- trimws(opts$output) 
if (is.null(opts$output)) {
    opts$output <- paste0("output/phylosig_pa_", tree_file_name, "_lam", opts$lam, "_db", opts$decision_boundary, "_draws", opts$N_draws, "_samples", opts$N_samples, ".rds")
}


saveRDS(res, file=opts$output)
