source("../reference.R")

### sufficiently small for null data
test_that("null data gives small Bayes factor", {
  for (seed in 2026:(2026 + 4)) {
    set.seed(seed)
    trees <- lapply(1:5, function(x) ape::rcoal(100))
    x <- setNames(rnorm(100), trees[[1]]$tip.label)
    bf <- lambdaBF(trees, x, N_samples = 10, return_trace = F, importance_sampling = T,
                   a = 3, d = 5)
    expect_lt(bf$bf, 1)
  }
})

### sufficiently large for synthetic data
test_that("synthetic data gives large Bayes factor", {
  for (seed in 2026:(2026 + 4)) {
    set.seed(seed)
    tree <- ape::rcoal(100)
    x <- phytools::fastBM(tree)
    bf <- lambdaBF(tree, x, N_samples = 100, return_trace = F, importance_sampling = T,
                   a = 3, d = 5)
    expect_gt(bf$bf, 1)
  }
})

### importance sampling gives similar results to non-importance sampling
test_that("importance sampling gives similar results to non-importance sampling", {
  for (seed in 2026:(2026 + 4)) {
    set.seed(seed)
    tree <- ape::rcoal(150)
    x <- phytools::fastBM(tree)
    bf_importance <- lambdaBF(tree, x, N_samples = 75, return_trace = F, importance_sampling = T,
                             a = 3, d = 5)
    bf_no_importance <- lambdaBF(tree, x, N_samples = 250, return_trace = F, importance_sampling = F,
                                a = 3, d = 5)
    expect_equal(log10(bf_importance$bf), log10(bf_no_importance$bf), tolerance=0.005)
  }
})


### compare with former implementation using slow likelihood
test_that("new implementation gives similar results to former implementation", {
  for (seed in 2026:(2026 + 4)) {
    set.seed(seed)
    tree <- ape::rcoal(100)
    x <- phytools::fastBM(tree)
    bf_new <- lambdaBF(tree, x, N_samples = 250, return_trace = F, importance_sampling = T,
                      a = 3, d = 5)
    bf_old <- lambdaBF_old(tree, x, N_samples = 250, return_trace = F, importance_sampling = T,
                          a = 3, d = 5)
    expect_equal(log10(bf_new$bf), log10(bf_old$bf), tolerance=0.005)
  }
})