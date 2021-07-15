make_knots <- function(M, n_knots, start, end) {
    
    if(!(M == length(n_knots))){
        stop('Specify number of knots for each resolution')
    }
    
    knots   <- vector(mode = "list", length = M)
    for(i in 1:M) {
        knots[[i]] <- seq(start, end, length.out = n_knots[i])
    }
    
    return(knots)
}

make_kernel <- function(x, M, knots, radius, form = "gaussian") {
    # evaluate the kernel/basis in 1d
    if (!(form %in% c("uniform", "gaussian", "exponential", "polynomial"))) {
        stop('form must be either "uniform" or "gaussian" or "exponential" or "polynomial" ...')
    }
    
    n <- length(x)
    kernel_list <- vector(mode = "list", length = M)
    
    for(i in 1:M) {
        n_knots <- length(knots[[i]])
        kernel_list[[i]] <- matrix(0, n, n_knots)
        d <- fields::rdist(x, knots[[i]])
        
        if (form == "uniform") {
            kernel_list[[i]] <- ifelse(d < radius[i], 1, 0)
        }
        if (form == "gaussian") {
            kernel_list[[i]] <- dnorm(d / radius[i])
        }
        if (form == "polynomial") {
            kernel_list[[i]] <- ifelse(d < radius[i], (1 - d^3 / radius[i]^3)^3, 0)
        }
        if (form == "exponential") {
            kernel_list[[i]] <- dexp(d / radius[i])
        }
    }
    
    kernel <- do.call(cbind, kernel_list)
    
    return(kernel)
}

mcmc <- function(y, x, 
                 M, knots, radius, form = "gaussian",
                 n_mcmc = 5000, burnin = 2500, n_message = 500) {
    
    N <- length(y)
    n_knots <- unlist(lapply(knots, length))
    total_n_knots <- sum(n_knots)
    
    ## get prior parameters
    a_epsilon0 <- 0.01
    b_epsilon0 <- 0.01
    a_alpha0 <- 0.01
    b_alpha0 <- 0.01
    
    ## initialize parameters
    sigma2_epsilon <- runif(1, 1, 5)
    sigma2_alpha <- rep(1, M)
    alpha <- rnorm(total_n_knots)
    Z <- make_kernel(x, M, knots, radius, form = "gaussian")
    
    ## precalculate values
    tZZ <- t(Z) %*% Z
    tZy <- t(Z) %*% y
    I_m <- diag(total_n_knots)
    
    dims_idx <- rep(1, n_knots[1])
    if(M > 1) {
        for(m in 2:M) {
            dims_idx <- c(dims_idx, rep(m, n_knots[m]))
        }
    }
    
    ## set up save variables
    sigma2_epsilon_save <- rep(0, n_mcmc-burnin)
    sigma2_alpha_save <- matrix(0, n_mcmc-burnin, M)
    alpha_save <- matrix(0, n_mcmc-burnin, total_n_knots)
    
    ## mcmc loop
    for(k in 1:n_mcmc) {
        # if (k%% n_message == 0) {
        #     message("Iteration ", k , " out of ", n_mcmc)
        # }
        
        ## sample sigma2 alpha
        for(m in 1:M) {
            alpha_m <- alpha[dims_idx == m]
            sigma2_alpha[m] <- 1 / rgamma(1, a_alpha0 + n_knots[m]/2,
                                          b_alpha0 + (t(alpha_m) %*% alpha_m) / 2)
        }
        
        ## sample sigma2 epsilon
        sigma2_epsilon <- 1 / rgamma(1, a_epsilon0 + N/2,
                                     b_epsilon0 + (t(y - Z %*% alpha) %*% 
                                                       (y - Z %*% alpha))/2)
        
        
        ## sample alpha
        sigma2_alpha_inv <- rep(1/sigma2_alpha[1], n_knots[1])
        if(M > 1) {
            for(m in 2:M) {
                sigma2_alpha_inv <- c(sigma2_alpha_inv, 
                                      rep(1/sigma2_alpha[m], n_knots[m]))
            }
        }
        
        A_inv <- solve(tZZ / sigma2_epsilon + sigma2_alpha_inv*I_m)
        b <- tZy  / sigma2_epsilon
        alpha <- c(mvnfast::rmvn(1, A_inv %*% b, A_inv))
        
        ## save MCMC variables
        if(k > burnin) {
            i <- k-burnin
            sigma2_epsilon_save[i] <- sigma2_epsilon
            sigma2_alpha_save[i, ] <- sigma2_alpha
            alpha_save[i, ] <- alpha
        }
    }
    
    return(list(sigma2_epsilon=sigma2_epsilon_save,
                sigma2_alpha=sigma2_alpha_save,
                alpha=alpha_save))
}