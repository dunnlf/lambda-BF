### Applying lambda-BF to example data from papers
###
### Data from the following are considered:
###    - Gônet, J. et al. (2023)
###    - Zhang, ZN., Lyu, X., Niu, SZ. et al. (2026)
###    - Round, E., Dockum, R., Ryder, R. J. (2022)
###


library(ape)
library(phytools)
library(lambdaBF)

# ---- Gônet et al. ----

# load data, duplicate S. punctatus and E. elgans samples removed
data <- read.csv("lambda-BF-example-data\\gonet-et-al\\femur-measurements.csv", header = TRUE, row.names = 1)
species <- rownames(data)

# read trees
trees <- read.tree("lambda-BF-example-data\\gonet-et-al\\Phylo_Tree_Set.tre")

# compute/load bayes factors
bf_file <- "lambda-BF-example-data\\gonet-et-al\\bf_res.rData"
if (file.exists(bf_file)) {
    load(bf_file)
} else {
    bf_res <- list()
    for (i in 1:ncol(data)) {
        print(paste("Computing BF for ", colnames(data)[i]))
        bf_res[[i]] <- lambdaBF(trees, setNames(data[,i],
                                sapply(species, function(x) gsub(" ", "_", x))))
    }

    save(bf_res, file = bf_file)
}

# print bayes factors
bfs <- sapply(bf_res, function(x) x$bf)
log10(bfs)



