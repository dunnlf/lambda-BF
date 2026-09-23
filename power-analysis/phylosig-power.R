### power analysis of mean p-value from phylosig

library(ape)
library(phytools)
library(lambdaBF)

# read trees
trees_50 <- read.nexus("lambda-BF-example-data\\power-analysis\\mammal-trees-50.nex")
trees_100 <- read.nexus("lambda-BF-example-data\\power-analysis\\mammal-trees-100.nex")
trees_150 <- read.nexus("lambda-BF-example-data\\power-analysis\\mammal-trees-150.nex")
trees_250 <- read.nexus("lambda-BF-example-data\\power-analysis\\mammal-trees-250.nex")


# function to generate synthetic data
get_synthetic_data <- function(lam, n_samps, t) {
  n_species <- Ntip(t)
  species_names <- t$tip.label
  x_m1 <- matrix(nrow=n_species, ncol=n_samps)

  for (i in 1:n_samps){

    x_m1[,i] <- rnorm(n_species, mean=0, sd=1)
    vcv <- get_pagel_cov(lam, t)
    eig_vcv <- eigen(vcv)
    Q <- eig_vcv$vectors
    eigvals <- diag(eig_vcv$values)
    S = Q %*% sqrt(eigvals) %*% t(Q)
    
    # transform to have desired covariance structure
    x_m1[,i] <- S %*% x_m1[,i]
  }
  x_m1 <- data.frame(x_m1)
  rownames(x_m1) <- species_names
  return(x_m1)
}


lam_range <- seq(0.1, 1, by=0.1)

n_trees=50
n_draws=50
t_50 <- trees_50[1:n_trees]
t_100 <- trees_100[1:n_trees]
t_150 <- trees_150[1:n_trees]
t_250 <- trees_250[1:n_trees]

trees_list <- list(t_50, t_100, t_150, t_250)
name_list <- c("50", "100", "150", "250")


for (k in 1:4){
    print(paste("Running for", name_list[k], "species"))
    t <- trees_list[[k]]
    name <- name_list[k]

    res_lam <- matrix(ncol=length(lam_range), nrow=n_draws)
    
    for (i in 1:length(lam_range)){
        print(paste("Running lambda =", lam_range[i]))

        for (j in 1:n_draws){
            mcc <- phangorn::maxCladeCred(t)
            trait_data <- as.numeric(get_synthetic_data(
                lam = lam_range[i],
                n_samps = 1,
                t = mcc
            )[, 1])

            names(trait_data) <- mcc$tip.label

            mean_p <- average_lambda_pvalue(
                trees = t,
                trait_data = trait_data
            )

            res_lam[j, i] <- mean_p
            res_decision[j, i] <- as.integer(mean_p < 0.05)
        }
    }

    filename <- paste0(
        "lambda-BF-example-data/power-analysis/res_lambda_pvalue_",
        name, "_", n_draws, ".RData"
    )
    save(res_decision, res_lam, file = filename)

}