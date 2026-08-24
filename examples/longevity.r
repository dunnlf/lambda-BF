### Applying lambda-BF to example data from papers
###
### Data from the following are considered:
###    - Venditti, C., Baker, J. & Barton, R.A. (2024)
###    - Zhang, ZN., Lyu, X., Niu, SZ. et al. (2026)
###    - Chen, Z., Wiens, J.J. (2020)
###


library(ape)
library(phytools)
library(lambdaBF)

# ---- Zhang et al. ----

# read lifespan data
data <- read.csv("lambda-BF-example-data\\zhang\\zhang-lifespan.csv")

species <- data$Species

# write species to csv
write.csv(species, "lambda-BF-example-data\\zhang\\zhang-species.csv", row.names = FALSE)

# read tree from timetree
tree <- read.tree("lambda-BF-example-data\\zhang\\zhang-species.tre")

# bayes factor on tree used by Zhang et al. (2026)
LQ <- setNames(data$Longevity.quotient..LQ., sapply(species,
                                function(x) gsub(" ", "_", x)))

if (length(setdiff(names(LQ), tree$tip.label)) != 0) {
    print("Warning: some species in LQ data not found in trees")}

bf <- lambdaBF(tree, LQ)


# on a sample of 100 trees from VertLife
trees_100 <- read.nexus("lambda-BF-example-data\\zhang\\mammal-trees-100.nex")

if (length(setdiff(names(LQ), trees_100[[1]]$tip.label)) != 0) {
    print("Warning: some species in LQ data not found in trees")}

bf_100 <- lambdaBF(trees_100, LQ, return_trace = TRUE)
