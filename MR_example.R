MR_model <- function() {
    library(ggplot2)
    library(gridExtra)
    library(latex2exp)
    source("make_kernel.R")
    source("mcmc_mr.R")
    
    ## set ggplot theme
    theme_set(theme_minimal())
    
    ## Simulate data ---------------------------------------------------------------
    N <- 50
    x <- seq(0, 1, length.out = N)
    n_knots1 <- 10
    n_knots2 <- 20
    n_knots3 <- 40
    
    knots1 <- seq(0, 1, length.out = n_knots1)
    knots2 <- seq(0, 1, length.out = n_knots2)
    knots3 <- seq(0, 1, length.out = n_knots3)
    knots_mr <- list(
        knots1 = knots1,
        knots2 = knots2,
        knots3 = knots3
    )
    
    sigma2_epsilon <- 0.5^2
    sigma2_alpha1 <- 1
    sigma2_alpha2 <- 1
    sigma2_alpha3 <- 1
    
    alpha1 <- rnorm(n_knots1, 0, sqrt(sigma2_alpha1))
    alpha2 <- rnorm(n_knots2, 0, sqrt(sigma2_alpha2))
    alpha3 <- rnorm(n_knots3, 0, sqrt(sigma2_alpha3))
    
    stddev <- c(1, 0.5^2, 0.25^2)
    Z1 <- make_kernel(x, knots1, stddev = sqrt(stddev[1]))
    Z2 <- make_kernel(x, knots2, stddev = sqrt(stddev[2]))
    Z3 <- make_kernel(x, knots3, stddev = sqrt(stddev[3]))
    
    z <- Z1 %*% alpha1 + Z2 %*% alpha2 + Z3 %*% alpha3
    y <- z + rnorm(N, 0, sqrt(sigma2_epsilon))
    
    plot_data <- data.frame(
        x = x,
        z = z,
        y = y
    )
    
    ggplot(data = plot_data) + 
        geom_point(aes(x=x, y=y)) + 
        geom_line(aes(x=x, y=z))
    
    ## MR model --------------------------------------------------------------------
    n_mcmc <- 10000
    burnin <- 5000
    mr_fit <- mcmc_mr(y, x, 
                      knots_mr, stddev=stddev, 
                      n_mcmc = n_mcmc, burnin = burnin, n_message = 500)
    
    
    
    ## Plot ------------------------------------------------------------------------
    plot_mcmc <- data.frame(
        Samples = 1:(n_mcmc-burnin),
        sigma2_epsilon = mr_fit$sigma2_epsilon,
        sigma2_alpha1 = mr_fit$sigma2_alpha1,
        sigma2_alpha2 = mr_fit$sigma2_alpha2,
        sigma2_alpha3 = mr_fit$sigma2_alpha3
    )
    
    p1 <- ggplot() + 
        geom_line(data = plot_mcmc, aes(x=Samples, y=sigma2_epsilon)) + 
        geom_hline(yintercept = sigma2_epsilon, color = "red") +
        ylab(TeX(r'($\sigma^2_{\epsilon}$)'))
    
    p2 <- ggplot() + 
        geom_line(data = plot_mcmc, aes(x=Samples, y=sigma2_alpha1)) + 
        geom_line(data = plot_mcmc, aes(x=Samples, y=sigma2_alpha2)) + 
        geom_line(data = plot_mcmc, aes(x=Samples, y=sigma2_alpha3)) + 
        geom_hline(yintercept = sigma2_alpha1) +
        geom_hline(yintercept = sigma2_alpha2) +
        geom_hline(yintercept = sigma2_alpha3) +
        ylab(TeX(r'($\sigma^2_{\alpha}$)'))
    
    plot_data <- data.frame(
        x = x,
        z = z,
        y = y,
        MR = Z1 %*% apply(mr_fit$alpha1, 2, mean) + 
            Z2 %*% apply(mr_fit$alpha2, 2, mean) + 
            Z3 %*% apply(mr_fit$alpha3, 2, mean)
    )
    
    p3 <- ggplot(data = plot_data) + 
        geom_point(aes(x=x, y=y)) + 
        geom_line(aes(x=x, y=z, color = "True")) + 
        geom_line(aes(x=x, y=MR, color = "MR")) + 
        theme(legend.title = element_blank()) + 
        theme(legend.position="bottom")
    
    grid.arrange(p1, p2, p3, nrow = 3)
}

