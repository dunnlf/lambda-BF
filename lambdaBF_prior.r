get_marginal_lhood <- function(lam, tau, x, logarithm=F, a, d) {
    # return factor-reduced marginal likelihood of Pagel's lambda model
    # facotrs independent of lambda, tau, x are removed
    # a, d are prior parameters of IG prior on variance
    # it is assumed that conditional prior of mean is N(0, a/(d-2))

    n <- length(x)
    C_lam <- get_pagel_cov(lam, tau)
    C_inv <- solve(C_lam)
    S <- sum(C_inv)

    x <- x[rownames(C_lam)]
    if (logarithm) {
        return(-0.5 * log(1+S) + (-0.5*(d+n))*log(a + t(x)%*%C_inv%*%x - ((sum(C_inv%*%x))**2)/(1+S)))
    }
    else {
    return(sqrt(1/(1+S)) * (a + t(x)%*%C_inv%*%x - ((sum(C_inv%*%x))**2)/(1+S))**(-0.5*(d+n)))
    }
}

get_full_marginal_lhood <- function(lam, tau, x, logarithm=F, a, d) {
    # return factor-reduced marginal likelihood of Pagel's lambda model
    # facotrs independent of lambda, tau, x are removed
    # a, d are prior parameters of IG prior on variance
    # it is assumed that conditional prior of mean is N(0, a/(d-2))

    n <- length(x)
    C_lam <- get_pagel_cov(lam, tau)
    C_inv <- solve(C_lam)
    S <- sum(C_inv)

    x <- x[rownames(C_lam)]

    if (logarithm) {
        return((0.5*d)*log(a) + log(gamma(0.5*(d+n))) -(n/2)*log(pi) - log(gamma(d)/2) -0.5 * log(1+S) +
        (-0.5*(d+n))*log(a + t(x)%*%C_inv%*%x - ((sum(C_inv%*%x))**2)/(1+S)))
    }
    else {
    return(a^(0.5*d) / pi^(n/2) * gamma(0.5*(d+n)) / gamma(d/2) * sqrt(1/(1+S)) * (a + t(x)%*%C_inv%*%x -
    ((sum(C_inv%*%x))**2)/(1+S))**(-0.5*(d+n)))
    }
}


lambdaBF_prior <- function (trees, x, N_samples = 100, return_trace = F, importance_sampling = T,
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

    for (t in trees) {
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

        if(is.ultrametric(t_i) == FALSE) {
            warning("Warning: Tree not ultrametric")
        }
        tree_age <- as.numeric(diag(vcv.phylo(t_i)))
        x_i <- sqrt(tree_age) * (x_i - mean(x_i)) / sd(x_i)

        log_lhood_func <- function(z) {
            get_marginal_lhood(z, t_i, x_i, logarithm = T, a = a, d = d)
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
        lhood_null_i <- get_marginal_lhood(0, t_i, x_i, logarithm = T, a = a, d = d)
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
