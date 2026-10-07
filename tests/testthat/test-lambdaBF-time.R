### timing lambdaBF
paste('Running timing...')
# time runs on 50, 500, 1000 tips, for 10, 50, and 100 samples
test_that("timing 50 tips", {
  set.seed(2026)
  tree <- ape::rcoal(50)
  cat("50 tips, 15 draws \n")
  for (n_samples in c(10, 50, 100)) {
    times_50 <- replicate(15, {
      x <- phytools::fastBM(tree)
      system.time(lambdaBF(tree, x, N_samples = n_samples, return_trace = F,
                           importance_sampling = T, a = 3, d = 5))["elapsed"]
    })

    cat(n_samples,
        "samples:", mean(times_50), "seconds\n")
  }
  succeed()
})

test_that("timing 250 tips", {
  set.seed(2026)
  tree <- ape::rcoal(250)
  cat("250 tips, 15 draws \n")
  for (n_samples in c(10, 50, 100)) {
    times_50 <- replicate(15, {
      x <- phytools::fastBM(tree)
      system.time(lambdaBF(tree, x, N_samples = n_samples, return_trace = F,
                           importance_sampling = T, a = 3, d = 5))["elapsed"]
    })

    cat(n_samples,
        "samples:", mean(times_50), "seconds\n")
  }
  succeed()
})


test_that("timing 500 tips", {
  set.seed(2026)
  tree <- ape::rcoal(500)
  cat("500 tips, 15 draws \n")
  for (n_samples in c(10, 50, 100)) {
    times_50 <- replicate(15, {
      x <- phytools::fastBM(tree)
      system.time(lambdaBF(tree, x, N_samples = n_samples, return_trace = F,
                           importance_sampling = T, a = 3, d = 5))["elapsed"]
    })

    cat(n_samples,
        "samples:", mean(times_50), "seconds\n")
  }
  succeed()
})



