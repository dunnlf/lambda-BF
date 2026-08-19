### function to evaluate evidence for lambda (Pagel, 1999) via the Bayes factor
### allowing for uncertainty in the phylogenetic tree
### written by Loki Dunn 2026


#' Approximate Bayes factor for Pagel model
#' 
#' Uses importance sampling to approximate the Bayes factor of Pagel's lambda
#' (Pagel, 1999), aggregating across trees. Importance sampling used to
#' marginalise over lambda and makes use of a Laplace approximation to the
#' likelihood.
#' For some trees, the likelihood may not have a maximum in (0,1), in which case 
#' the Laplace approximation fails. In this case lambda is sampled uniformly.
#'  
#' @param trees phylo object, or list of phylo objects
#' @param x named vector of observations for each tip
#' @param N_samples number of importance samples to draw
#' @param logarithm boolean, whether to return log likelihood.
#' @param return_trace boolean, whether to return sampling traces.
#' @param importance_sampling boolean, whether to use importance sampling.
#'
#' @return list containing Bayes factor estimate 'bf', as well as bayes factors
#' for each tree in 'bf_samples', optionally contains sampling traces
#' 
#' 
#' @examples
#' trees <- ape::rmtree(5, 100)
#' x <- phytools::fastBM(trees[[1]])
#' pagelBF(trees, x)
#' 
#' 
#' @export
pagelBF <- function(trees, x, N_samples=100, return_trace=F,
                                importance_sampling=T) {
  
  if (class(trees) == 'phylo'){  # if single tree provided wrap in list
    trees <- list(tree=trees)
  }
  
  if (is.null(names(x))){
    warning('Tip data not named, assuming same order as tree tips')
    names(x) <- trees[[1]]$tip.label
  }

  # species mismatch debugging
  num_dropped_species <- c()
  dropped_species <- list()
  
  # samples per tree
  pagel_samples <- c()
  null_samples <- c()
  bf_samples <- c()
  
  # traces
  trace <- c()
  pagel_trace <- list()
  null_trace <- list()
  
  # sampling methods, used for debugging Laplace approximate failures
  sampling_methods <- c()
  
  for (t in trees) {
    # check name mismatches and prune if necessary
    check_t <- geiger::name.check(t, x, data.names = names(x))
    
    if (length(check_t) != 1){
      t_i <- drop.tip(t, check_t$tree_not_data)
      x_i <- x[setdiff(names(x), check_t$data_not_tree)]
      num_dropped_species <- c(num_dropped_species,
                               length(check_t$data_not_tree) + length(check_t$tree_not_data))
      dropped_species[[length(dropped_species)+1]] <- c(check_t$data_not_tree,
                                                        check_t$tree_not_data)
    } else {
      x_i <- x
      t_i <- t
      num_dropped_species <- c(num_dropped_species, 0)
      dropped_species[[length(dropped_species)+1]] <- NULL
    }
  
    lhood_func <- function(z) {get_pagel_lhood(z, t_i, x_i)}
    
    # compute Laplace approximation to likelihood
    log_lhood_func <- function(z) {log(lhood_func(z))}
    
    lam_0 <- optimise(lhood_func, c(0,1), maximum=T)$maximum
    
    logL_2nd_div <- abs(numDeriv::hessian(log_lhood_func, lam_0))
    
    lap_approx <- function(z) {dnorm(z, mean=lam_0, sd=1/sqrt(logL_2nd_div))}
    
    # sampling
    if (is.nan(logL_2nd_div) | lam_0 > 0.99 | lam_0 < 0.01 | (!importance_sampling)) {
      # no maximum in (0,1), use uniform sampling
      lam_samps <- runif(N_samples)
      log_num <- c()
      log_den <- c()
      for (lam in lam_samps){
        
        # log-sum-exp trick to prevent underflow
        log_num <- c(log_num, log_lhood_func(lam))
        log_den <- c(log_den, 0) # uniform samples used
      }
      sampling_methods <- c(sampling_methods, 'unif')
      
    } else {
      # importance sampling
      lam_samps <- truncnorm::rtruncnorm(N_samples, a=0, b=1, mean=lam_0,
                              sd=1/sqrt(logL_2nd_div))
      
      log_num <- c()
      log_den <- c()
      for (lam in lam_samps){
        # log-sum-exp trick to prevent underflow
        log_num <- c(log_num, log(lhood_func(lam)))
        log_den <- c(log_den, log(lap_approx(lam)))
      }
      
      sampling_methods <- c(sampling_methods, 'imp')
    }
    
    ponent = max(log_num - log_den)
    exp_terms = exp(log_num - log_den - ponent)
    
    lhood_pagel_i <- exp(ponent) * mean(exp_terms)
    
    # null likelihood - no dependence on lambda
    lhood_null_i <- get_pagel_lhood(0, t_i, x_i)
    
    pagel_samples <- c(pagel_samples, lhood_pagel_i)
    null_samples <- c(null_samples, lhood_null_i)
    
    bf_samples <- c(bf_samples, lhood_pagel_i/lhood_null_i)
    
    if (return_trace) {
      trace <- c(trace,
                   log_sum_exp_mean(pagel_samples)/log_sum_exp_mean(null_samples))
      null_trace[[length(null_trace)+1]] <- get_trace(null_samples)
      pagel_trace[[length(pagel_trace)+1]] <- get_trace(exp(ponent) * exp_terms)
    }
  }
  
  # compute Bayes factor
  num <- log_sum_exp_mean(pagel_samples)
  den <- log_sum_exp_mean(null_samples)
  bf <- num/den
    
  
  if (return_trace){
    return(list("bf"=bf,
                "bf_samples"=bf_samples,
                "trace"=trace,
                "pagel_trace"=pagel_trace,
                "null_trace"=null_trace,
                "samplers"=sampling_methods,
                "num_dropped_species"=num_dropped_species,
                "dropped_species"=dropped_species,
                "pagel_mls"=pagel_samples,
                "null_mls"=null_samples))
  } else {
    return(list('bf'=bf,
                'bf_samples'=bf_samples))
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


#' Compute pagel likelihood
#' 
#' For a given value of lambda, computes the likelihood of the Pagel model
#'  
#' @param lam scalar value.
#' @param tau tree object of class 'phylo'.
#' @param x trait data for tips of tau.
#' @param logarithm boolean, whether to return log likelihood.
#'
#' @return scalar value, likelihood of pagel model given lam
#' 
#' 
#' @examples
#' t <- ape::rtree(10)
#' x <- phytools::fastBM(t)
#' get_pagel_lhood(1, t, x)
#' 
#' 
#' @export
get_pagel_lhood <- function(lam, tau, x, logarithm=F) {
  
  C_lam_tau <- get_pagel_cov(lam, tau)
  C_inv <- solve(C_lam_tau)
  n <- nrow(C_inv)
  
  # align data with matrix columns
  x_ <- x[rownames(C_lam_tau)]
  x_ <- matrix(x_, ncol=1)
  
  # estimate mean and variance
  mu <- as.numeric(sum(C_inv %*% x_)/sum(C_inv))
  sig2 <- as.numeric(t(x_ - mu) %*% C_inv %*% (x_ - mu)/n)
  
  if(logarithm){
    (-n/2) * log(2 * pi) - (1/2) * determinant(sig2*C_lam_tau, logarithm = T)$modulus - 
      (1/2) * t(x_ - mu) %*% (C_inv/sig2) %*% (x_ - mu)
    
  } else {
    (2*pi)**(-n/2) * det(sig2 * C_lam_tau)**(-1/2) *
      exp(-(1/2) * t(x_ - mu) %*% (C_inv/sig2) %*% (x_ - mu))
    
  }
}
