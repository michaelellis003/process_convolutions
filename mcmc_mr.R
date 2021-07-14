mcmc_mr <- function(y, x, 
                    knots, stddev, 
                    n_mcmc = 5000, burnin = 2500, n_message = 500) {
    
    N <- length(y)
    m1 <- length(knots$knots1)
    m2 <- length(knots$knots2)
    m3 <- length(knots$knots3)
    m <- m1 + m2 + m3
    
    ## get prior parameters
    a_epsilon0 <- 0.01
    b_epsilon0 <- 0.01
    a_alpha0 <- 0.01
    b_alpha0 <- 0.01
    
    ## initialize parameters
    sigma2_epsilon1 <- runif(1, 1, 5)
    
    sigma2_alpha1 <- runif(1, 1, 5)
    sigma2_alpha2 <- runif(1, 1, 5)
    sigma2_alpha3 <- runif(1, 1, 5)
    
    alpha1 <- rnorm(m1, 0, 1)
    alpha2 <- rnorm(m2, 0, 1)
    alpha3 <- rnorm(m3, 0, 1)
    
    Z1 <- make_kernel(x, knots$knots1, stddev[1])
    Z2 <- make_kernel(x, knots$knots2, stddev[2])
    Z3 <- make_kernel(x, knots$knots3, stddev[3])
    
    ## precalculate values
    tZZ1 <- t(Z1) %*% Z1
    tZZ2 <- t(Z2) %*% Z2
    tZZ3 <- t(Z3) %*% Z3
    
    tZ1y <- t(Z1) %*% y
    tZ2y <- t(Z2) %*% y
    tZ3y <- t(Z3) %*% y
    
    I_m1 <- diag(m1)
    I_m2 <- diag(m2)
    I_m3 <- diag(m3)
    
    ## set up save variables
    sigma2_epsilon_save <- rep(0, n_mcmc-burnin)
    
    sigma2_alpha1_save <- rep(0, n_mcmc-burnin)
    sigma2_alpha2_save <- rep(0, n_mcmc-burnin)
    sigma2_alpha3_save <- rep(0, n_mcmc-burnin)
    
    alpha1_save <- matrix(0, n_mcmc-burnin, m1)
    alpha2_save <- matrix(0, n_mcmc-burnin, m2)
    alpha3_save <- matrix(0, n_mcmc-burnin, m3)
    
    ## mcmc loop
    for(k in 1:n_mcmc) {
        # if (k%% n_message == 0) {
        #     message("Iteration ", k , " out of ", n_mcmc)
        # }
        
        ## sample sigma2 alpha 1
        sigma2_alpha1 <- 1 / rgamma(1, a_alpha0 + m1/2,
                                    b_alpha0 + (t(alpha1) %*% alpha1) / 2)

        ## sample sigma2 alpha 2
        sigma2_alpha2 <- 1 / rgamma(1, a_alpha0 + m2/2,
                                    b_alpha0 + (t(alpha2) %*% alpha2) / 2)

        ## sample sigma2 alpha 3
        sigma2_alpha3 <- 1 / rgamma(1, a_alpha0 + m3/2,
                                    b_alpha0 + (t(alpha3) %*% alpha3) / 2)
        
        ## sample sigma2 epsilon
        Zalpha <- Z1 %*% alpha1 + Z2 %*% alpha2 + Z3 %*% alpha3
        sigma2_epsilon <- 1 / rgamma(1, a_epsilon0 + N/2,
                                     b_epsilon0 + (t(y - Zalpha) %*% (y - Zalpha))/2)
        
        
        ## sample alpha1
        A_inv1 <- solve(tZZ1 / sigma2_epsilon + (1/sigma2_alpha1)*I_m1)
        b1 <- tZ1y  / sigma2_epsilon
        alpha1 <- c(mvnfast::rmvn(1, A_inv1 %*% b1, A_inv1))
        
        ## sample alpha2
        A_inv2 <- solve(tZZ2 / sigma2_epsilon + (1/sigma2_alpha2)*I_m2)
        b2 <- tZ2y  / sigma2_epsilon
        alpha2 <- c(mvnfast::rmvn(1, A_inv2 %*% b2, A_inv2))
        
        ## sample alpha3
        A_inv3 <- solve(tZZ3 / sigma2_epsilon + (1/sigma2_alpha3)*I_m3)
        b3 <- tZ3y  / sigma2_epsilon
        alpha3 <- c(mvnfast::rmvn(1, A_inv3 %*% b3, A_inv3))
        
        ## save MCMC variables
        if(k > burnin) {
            i <- k-burnin
            sigma2_epsilon_save[i] <- sigma2_epsilon
            sigma2_alpha1_save[i] <- sigma2_alpha1
            sigma2_alpha2_save[i] <- sigma2_alpha2
            sigma2_alpha3_save[i] <- sigma2_alpha3
            alpha1_save[i, ] <- alpha1
            alpha2_save[i, ] <- alpha2
            alpha3_save[i, ] <- alpha3
        }
    }
    
    return(list(sigma2_epsilon = sigma2_epsilon_save,
                sigma2_alpha1 = sigma2_alpha1_save,
                sigma2_alpha2 = sigma2_alpha2_save,
                sigma2_alpha3 = sigma2_alpha3_save,
                alpha1 = alpha1_save,
                alpha2 = alpha2_save,
                alpha3 = alpha3_save))
}