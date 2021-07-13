mcmc_basic <- function(y, x, 
                 knots, stddev, 
                 n_mcmc = 5000, burnin = 2500, n_message = 500) {
    
    N <- length(y)
    m <- length(knots)
    
    ## get prior parameters
    a_epsilon0 <- 0.01
    b_epsilon0 <- 0.01
    a_alpha0 <- 0.01
    b_alpha0 <- 0.01
    
    ## initialize parameters
    sigma2_epsilon <- runif(1, 1, 5)
    sigma2_alpha <- runif(1, 1, 5)
    alpha <- rnorm(m)
    Z <- make_kernel(x, knots, stddev)
    
    ## precalculate values
    tZZ <- t(Z) %*% Z
    tZy <- t(Z) %*% y
    I_m <- diag(m)
    
    ## set up save variables
    sigma2_epsilon_save <- rep(0, n_mcmc-burnin)
    sigma2_alpha_save <- rep(0, n_mcmc-burnin)
    alpha_save <- matrix(0, n_mcmc-burnin, m)
    
    ## mcmc loop
    for(k in 1:n_mcmc) {
        # if (k%% n_message == 0) {
        #     message("Iteration ", k , " out of ", n_mcmc)
        # }
        
        ## sample sigma2 alpha
        sigma2_alpha <- 1 / rgamma(1, a_alpha0 + m/2,
                                   b_alpha0 + (t(alpha) %*% alpha) / 2)
        
        ## sample sigma2 epsilon
        sigma2_epsilon <- 1 / rgamma(1, a_epsilon0 + N/2,
                                     b_epsilon0 + (t(y - Z %*% alpha) %*% 
                                                       (y - Z %*% alpha))/2)
        
        
        ## sample alpha
        A_inv <- solve(tZZ / sigma2_epsilon + (1/sigma2_alpha)*I_m)
        b <- tZy  / sigma2_epsilon
        alpha <- c(mvnfast::rmvn(1, A_inv %*% b, A_inv))
        
        ## save MCMC variables
        if(k > burnin) {
            i <- k-burnin
            sigma2_epsilon_save[i] <- sigma2_epsilon
            sigma2_alpha_save[i] <- sigma2_alpha
            alpha_save[i, ] <- alpha
        }
    }
    
    return(list(sigma2_epsilon=sigma2_epsilon_save,
                sigma2_alpha=sigma2_alpha_save,
                alpha=alpha_save))
}