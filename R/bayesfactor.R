### function to evaluate evidence for lambda (Pagel, 1999) via the Bayes factor
### allowing for uncertainty in the phylogenetic tree
### written by Loki Dunn 2026


#' Approximate Bayes factor for Pagel model
#' 
#' Uses importance sampling to approximate the Bayes factor of Pagel's lambda
#' (Pagel, 1999), aggregating across trees. Importance sampling used to
#' marginalise over lambda and makes use of a Laplace approximation to the
#' likelihood. A normal-inverse-gamma prior is used for the mean and variance of the trait.
#' 
#' For some trees, the likelihood may not have a maximum in (0,1), in which case 
#' the Laplace approximation fails. In this case lambda is sampled uniformly.
#'  
#' @param trees phylo object, or list of phylo objects
#' @param x named vector of observations for each tip
#' @param N_samples number of importance samples to draw
#' @param return_trace boolean, whether to return sampling traces.
#' @param importance_sampling boolean, whether to use importance sampling.
#' @param a scale parameter for the variance prior.
#' @param d shape parameter for the variance prior.
#'
#' @return list containing Bayes factor estimate 'bf', as well as bayes factors
#' for each tree in 'bf_samples', optionally contains sampling traces
#' 
#' 
#' @examples
#' trees <- ape::rmtree(5, 100)
#' x <- phytools::fastBM(trees[[1]])
#' lambdaBF(trees, x)
#' 
#' 
#' @export
lambdaBF <-  function (trees, x, N_samples = 100, return_trace = F, importance_sampling = T,
                        a = 3, d = 5) {
    if (class(trees) == "phylo") {
        trees <- list(tree = trees)
    }
    if (is.null(names(x))) {
        warning("Tip data not named, assuming same order as tree tips")
        names(x) <- trees[[1]]$tip.label
    }
    num_dropped_species <- c()
    dropped_species <- list()
    pagel_samples <- c()
    null_samples <- c()
    bf_samples <- c()
    trace <- c()
    pagel_trace <- list()
    null_trace <- list()
    sampling_methods <- c()

    # for log-sum-exp trick
    exponent_list <- c()
    pagel_logs <- list()
    null_logs <- c()

    ultrametric_warn_flag <- 0

    for (i in seq_along(trees)) {
        t <- trees[[i]]
        check_t <- geiger::name.check(t, x, data.names = names(x))
        if (length(check_t) != 1) {
            t_i <- drop.tip(t, check_t$tree_not_data)
            x_i <- x[setdiff(names(x), check_t$data_not_tree)]

            num_dropped_species <- c(num_dropped_species, length(check_t$data_not_tree) + 
                length(check_t$tree_not_data))
            dropped_species[[length(dropped_species) + 1]] <- c(check_t$data_not_tree, 
                check_t$tree_not_data)
        }
        else {
            x_i <- x
            t_i <- t
            num_dropped_species <- c(num_dropped_species, 0)
            dropped_species[[length(dropped_species) + 1]] <- NULL
        }

        tree_age <- as.numeric(diag(vcv.phylo(t_i)))

        if(is.ultrametric(t_i, tol = 1e-6*min(tree_age)) == FALSE & ultrametric_warn_flag == 0){
            warning("Warning: Tree not ultrametric")
            ultrametric_warn_flag <- 1

        }

        x_i <- sqrt(tree_age) * (x_i - mean(x_i)) / sd(x_i)

        log_lhood_func <- function(z) {
            get_pagel_lhood(z, t_i, x_i, logarithm = T, a = a, d = d)
        }

        lam_0 <- optimise(log_lhood_func, c(0, 1), maximum = T)$maximum
        logL_2nd_div <- abs(numDeriv::hessian(log_lhood_func, 
            lam_0))
        lap_approx <- function(z) {
            dnorm(z, mean = lam_0, sd = 1/sqrt(logL_2nd_div))
        }
        if (is.nan(logL_2nd_div) | lam_0 > 0.99 | lam_0 < 0.01 | 
            (!importance_sampling)) {
            lam_samps <- runif(N_samples)
            log_num <- c()
            log_den <- c()
            for (lam in lam_samps) {
                log_num <- c(log_num, log_lhood_func(lam))
                log_den <- c(log_den, 0)
            }
            sampling_methods <- c(sampling_methods, "unif")
        }
        else {
            lam_samps <- truncnorm::rtruncnorm(N_samples, a = 0, 
                b = 1, mean = lam_0, sd = 1/sqrt(logL_2nd_div))
            log_num <- c()
            log_den <- c()
            for (lam in lam_samps) {
                log_num <- c(log_num, log_lhood_func(lam))
                log_den <- c(log_den, log(lap_approx(lam)))
            }
            sampling_methods <- c(sampling_methods, "imp")
        }
        # apply same exponent trick to numerator and denominator
        # store exponents and logs to apply log-sum-exp to final result also
        lhood_null_i <- get_pagel_lhood(0, t_i, x_i, logarithm = T, a = a, d = d)
        ponent = max(c(log_num - log_den, lhood_null_i))

        exp_terms = exp(log_num - log_den - ponent)
        exp_null = exp(lhood_null_i - ponent)

        pagel_logs[[length(pagel_logs) + 1]] <- log_num - log_den
        null_logs <- c(null_logs, lhood_null_i)

        lhood_pagel_i <- exp(ponent) * mean(exp_terms)
    
        exponent_list <- c(exponent_list, ponent)

        pagel_samples <- c(pagel_samples, log(lhood_pagel_i))
        null_samples <- c(null_samples, lhood_null_i)
        bf_samples <- c(bf_samples, mean(exp_terms)/exp_null)
        if (return_trace) {
            ponent <- max(exponent_list)

            num <- mean(sapply(pagel_logs, function(x) mean(exp(x - ponent))))
            den <- mean(exp(null_logs - ponent))
            trace <- c(trace, num/den)

            # traces within individual trees, for debugging
            null_trace[[length(null_trace) + 1]] <- get_trace(null_samples)
            pagel_trace[[length(pagel_trace) + 1]] <- get_trace(exp(ponent) * 
                exp_terms)
        }

        # warn of excessive variance in importance sampling, measured relative to mean on log scale
        if (sd(log_num - log_den)/sqrt(N_samples) > 0.1*abs(mean(log_num - log_den))) {
            warning(paste("Warning: High variance for tree", 
                i, " consider increasing importance sampling draws"))
        }

    }
    ponent <- max(exponent_list)

    num <- mean(sapply(pagel_logs, function(x) mean(exp(x - ponent))))
    den <- mean(exp(null_logs - ponent))

    bf <- num/den

    if (return_trace) {
        return(list(bf = bf, bf_samples = bf_samples, trace = trace, 
            pagel_trace = pagel_trace, null_trace = null_trace, 
            samplers = sampling_methods, num_dropped_species = num_dropped_species, 
            dropped_species = dropped_species, pagel_mls = pagel_samples, 
            null_mls = null_samples))
    }
    else {
        return(list(bf = bf, bf_samples = bf_samples))
    }
}



#' Compute trace of mean estimate
#'
#' For a vector of values, x, returns a vector whose nth entry is the mean of
#' the first n elements of x.
#' 
#' @param x numeric vector.
#'
#' @return numeric vector giving element-wise means along x
#' @noRd
get_trace <- function(x) {
  # 
  trace = c(x[1])
  i=1
  for (x_i in x[2:length(x)]) {
    
    trace <- c(trace, trace[length(trace)]*i/(i+1) + x_i/(i+1))
    i <- i+1
  }
  return(trace)
}


#' Log-sum-exp mean
#' 
#' Computes mean of a vector of values, using log-sum-exp trick
#' to prevent underflow
#'  
#' @param x numeric vector.
#'
#' @return scalar value, mean of x
#' @noRd
log_sum_exp_mean <- function(x){
  
  logs <- log(x)
  ponent <- max(log(x))
  norm_sum <- exp(logs - ponent)
  
  exp(ponent) * mean(norm_sum)
}


#' Get Pagel covariance
#' 
#' For a tree and value of Pagel's lambda (Pagel, 1999), returns brownian motion
#' covariance matrix with off-diagonal entries scaled by lambda.
#' 
#' @param lam scalar value.
#' @param tau tree object of class 'phylo'.
#' 
#' @return named matrix of size (Ntip(tau), Ntip(tau))
#' 
#' 
#' @examples
#' t <- ape::rtree(10)
#' get_pagel_cov(0.5, t)
#' 
#' 
#' @export
get_pagel_cov <- function(lam, tau) {
  
  sig_vcv <- ape::vcv.phylo(tau)
  diag_vcv <- diag(sig_vcv)
  
  sig_lam <- sig_vcv - diag(diag_vcv)
  out_vcv <- diag(diag_vcv) + lam*sig_lam
  
  out_vcv
}


#' Compute marginal Pagel likelihood
#' 
#' For a given value of lambda, computes the (factor-reduced) marginal likelihood of the Pagel model
#' when marginalising over a normal-inverse-gamma prior on the mean and variance
#'  
#' @param lam scalar value.
#' @param tau tree object of class 'phylo'.
#' @param x trait data for tips of tau.
#' @param a scale parameter for the variance prior.
#' @param d shape parameter for the variance prior.
#' @param logarithm boolean, whether to return log likelihood.
#'
#' @return scalar value, likelihood of pagel model given lam
#' 
#' @examples
#' t <- ape::rtree(10)
#' x <- phytools::fastBM(t)
#' get_pagel_lhood(1, t, x)
#' 
#' 
#' @export
get_pagel_lhood <- function(lam, tau, x, logarithm=F, a, d) {
    # return factor-reduced marginal likelihood of Pagel's lambda model
    # facotrs independent of lambda, tau, x are removed
    # a, d are prior parameters of IG prior on variance
    # it is assumed that conditional prior of mean is N(0, a/(d-2))

    n <- length(x)
    C_lam <- get_pagel_cov(lam, tau)
    C_inv <- solve(C_lam)
    S <- sum(C_inv)

    x <- x[rownames(C_lam)]

    C_inv_x <- C_inv %*% x

    if (logarithm) {
        return(-0.5 * log(1+S) + (-0.5*(d+n))*log(a + t(x)%*%C_inv_x - ((sum(C_inv_x))**2)/(1+S)))
    }
    else {
    return(sqrt(1/(1+S)) * (a + t(x)%*%C_inv_x - ((sum(C_inv_x))**2)/(1+S))**(-0.5*(d+n)))
    }
}
