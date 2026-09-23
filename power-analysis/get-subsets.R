### selecting subsets of taxa for power analysis

library(ape)
library(phytools)

set.seed(2026)

# read 100 trees, VertLife mammal trees for 927 species from Zhang, ZN., Lyu, X., Niu, SZ. et al. (2026)
trees <- read.nexus("lambda-BF-example-data\\power-analysis\\mammal_tree_samples.nex")

    
# select subsets of species
species_50 <- sample(trees[[1]]$tip.label, 50)
species_100 <- sample(trees[[1]]$tip.label, 100)
species_150 <- sample(trees[[1]]$tip.label, 150)
species_250 <- sample(trees[[1]]$tip.label, 250)

trees_50 <- lapply(trees, function(x) keep.tip(x, species_50))
trees_100 <- lapply(trees, function(x) keep.tip(x, species_100))
trees_150 <- lapply(trees, function(x) keep.tip(x, species_150))
trees_250 <- lapply(trees, function(x) keep.tip(x, species_250))

# Check that every subset tree is ultrametric and bifurcating.
tree_checks <- lapply(
	list(trees_50, trees_100, trees_150, trees_250),
	function(tree_set) {
		data.frame(
			ultrametric = vapply(tree_set, is.ultrametric, logical(1), tol=1e-6),
			bifurcating = vapply(tree_set, is.binary, logical(1))
		)
	}
)

names(tree_checks) <- c("trees_50", "trees_100", "trees_150", "trees_250")

stopifnot(vapply(tree_checks, function(x) all(x$ultrametric), logical(1)))
stopifnot(vapply(tree_checks, function(x) all(x$bifurcating), logical(1)))

write.nexus(trees_50, file = "lambda-BF-example-data\\power-analysis\\mammal-trees-50.nex")
write.nexus(trees_100, file = "lambda-BF-example-data\\power-analysis\\mammal-trees-100.nex")
write.nexus(trees_150, file = "lambda-BF-example-data\\power-analysis\\mammal-trees-150.nex")
write.nexus(trees_250, file = "lambda-BF-example-data\\power-analysis\\mammal-trees-250.nex")

write.csv(species_50, file = "lambda-BF-example-data\\power-analysis\\mammal-species-50.csv", row.names = FALSE)
write.csv(species_100, file = "lambda-BF-example-data\\power-analysis\\mammal-species-100.csv", row.names = FALSE)
write.csv(species_150, file = "lambda-BF-example-data\\power-analysis\\mammal-species-150.csv", row.names = FALSE)
write.csv(species_250, file = "lambda-BF-example-data\\power-analysis\\mammal-species-250.csv", row.names = FALSE)
