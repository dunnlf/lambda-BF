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
data <- subset(data, select = -c(n_forms, FL_VC))

bf_res_file <- "lambda-BF-example-data/round-et-al/bf_res.rData"
if (file.exists(bf_res_file)) {
        load(bf_res_file)
} else {
        bf_res <- list()
        for (i in 1:ncol(data)) {
                print(paste("Computing BF for ", colnames(data)[i]))
                bf_res[[i]] <- lambdaBF(trees, setNames(data[,i],
                                                                sapply(langs, function(x) gsub(" ", "_", x))))
        }

        save(bf_res, file = bf_res_file)
}

# print results to compare with those of paper
print(log10(bf_res[[1]]$bf))
hist(log10(bf_res[[1]]$bf_samples),
        main = 'FL_Vlength')

print(log10(bf_res[[3]]$bf))
hist(log(bf_res[[3]]$bf_samples),
        main ='FL_Cmanner')

print(log10(bf_res[[4]]$bf))
hist(log10(bf_res[[4]]$bf_samples),
        main ='FL_Cplace')

