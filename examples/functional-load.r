### Applying lambda-BF to example data from papers
###
### Data from the following are considered:
###    - Gônet, J. et al. (2023)
###    - Round, E., Dockum, R., Ryder, R. J. (2022)
###    - Arato, J., Fitch, W. T. (2021)
###


library(ape)
library(phytools)
library(lambdaBF)

# ---- Round et al. ----

# load functional load (FL) data
data <- read.csv("lambda-BF-example-data/round-et-al/FL_VlengthC_normalised.tsv", header = TRUE, row.names = 11,
            sep='\t')
langs <- rownames(data)

trees <- read.nexus("lambda-BF-example-data/round-et-al/PNY10_285.copy.trees")

# prune trees to data
trees <- lapply(trees, function(tree) {
    drop.tip(tree, setdiff(tree$tip.label, rownames(data)))
})
            
# apply lambda-BF to each FL measure in data
data <- subset(data, select = c(FL_Vlength, FL_Cmanner, FL_Cplace))

bf_res_file <- "lambda-BF-example-data/round-et-al/bf_res.rData"
if (file.exists(bf_res_file)) {
        load(bf_res_file)
} else {
        bf_res <- list()
        for (i in 1:ncol(data)) {
                print(paste("Computing BF for ", colnames(data)[i]))
                bf_res[[i]] <- lambdaBF(trees, setNames(data[,i], langs))
        }

        save(bf_res, file = bf_res_file)
}

# print results to compare with those of paper
print(paste("log10BF for FL_Vlength:", log10(bf_res[[1]]$bf)))
hist(log10(bf_res[[1]]$bf_samples),
        main = 'FL_Vlength')

print(paste("log10BF for FL_Cmanner:", log10(bf_res[[2]]$bf)))
hist(log10(bf_res[[2]]$bf_samples),
        main ='FL_Cmanner')

print(paste("log10BF for FL_Cplace:", log10(bf_res[[3]]$bf)))
hist(log10(bf_res[[3]]$bf_samples),
        main ='FL_Cplace')
