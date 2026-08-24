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

LQ <- setNames(data$Longevity.quotient..LQ., sapply(species,
                                function(x) gsub(" ", "_", x)))
bf <- lambdaBF(tree, LQ)
