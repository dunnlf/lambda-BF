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

# ---- Zhang et al. ----

# read lifespan data
data <- read.csv("lambda-BF-example-data\\zhang-et-al\\zhang-lifespan.csv")

species <- data$Species

# write species to csv
write.csv(species, "lambda-BF-example-data\\zhang-et-al\\zhang-species.csv", row.names = FALSE)

# read tree from timetree, used in paper
tree <- read.tree("lambda-BF-example-data\\zhang-et-al\\zhang-species.tre")

# bayes factor on tree used by Zhang et al. (2026)
LQ <- setNames(data$Longevity.quotient..LQ., sapply(species,
                                function(x) gsub(" ", "_", x)))

if (length(setdiff(names(LQ), tree$tip.label)) != 0) {
    print("Warning: some species in LQ data not found in trees")}

bf_file <- "lambda-BF-example-data\\zhang-et-al\\bf_res.rData"
if (file.exists(bf_file)) {
    load(bf_file)
} else {
    bf_res <- lambdaBF(tree, LQ)
    save(bf_res, file = bf_file)
}

# on a sample of 100 trees from VertLife
trees_100 <- read.nexus("lambda-BF-example-data\\zhang-et-al\\mammal-trees-100.nex")

if (length(setdiff(names(LQ), trees_100[[1]]$tip.label)) != 0) {
    print("Warning: some species in LQ data not found in trees")}

bf_file_100 <- "lambda-BF-example-data\\zhang-et-al\\bf_res_100.rData"
if (file.exists(bf_file_100)) {
    load(bf_file_100)
} else {
    bf_100 <- lambdaBF(trees_100, LQ)
    save(bf_100, file = bf_file_100)
}


# print summaries
