run_example <- function() {
    library(ggplot2)
    source("functions.R")
    
    theme_set(theme_minimal())
    
    ## Simulate data ---------------------------------------------------------------
    N <- 30
    x <- seq(0, 10, length.out = N)
    sigma2_epsilon <- 0.1^2
    Z_alpha <- sin(2*pi*(x/10)) + 0.2*sin(2*pi*(x/10))
    y <- Z_alpha + rnorm(N, 0, sqrt(sigma2_epsilon))
    
    ## Basic 1 resolution ----------------------------------------------------------
    M <- 1 # number of resolutions
    n_knots <- 7
    radius <- 2
    knots <- make_knots(M, n_knots, start = -1, end = 12)
    Z_basic <- make_kernel(x, M, knots, radius, form = "gaussian")
    
    n_mcmc <- 5000
    burnin <- 2500
    fit_basic <- mcmc(y, x, 
                      M, knots, radius = radius, form = "gaussian",
                      n_mcmc = n_mcmc, burnin = burnin, n_message = 500)
    
    
    ## MR model --------------------------------------------------------------------
    M <- 3 # number of resolutions
    n_knots <- c(7, 14, 28)
    radius <- c(2, 1, 0.5)
    knots <- make_knots(M, n_knots, start = -1, end = 12)
    Z_mr <- make_kernel(x, M, knots, radius, form = "gaussian")
    
    n_mcmc <- 5000
    burnin <- 2500
    fit_mr <- mcmc(y, x, 
                   M, knots, radius = radius, form = "gaussian",
                   n_mcmc = n_mcmc, burnin = burnin, n_message = 500)
    
    ## Plot ------------------------------------------------------------------------
    plot_data <- data.frame(
        x = x,
        z = Z_alpha,
        y = y,
        Basic = Z_basic %*% apply(fit_basic$alpha, 2, mean),
        MR = Z_mr %*% apply(fit_mr$alpha, 2, mean) 
    )
    
    ggplot(data = plot_data) + 
        geom_point(aes(x=x, y=y), size=2.5) + 
        geom_line(aes(x=x, y=z, color = "True"), lwd=1.5) + 
        geom_line(aes(x=x, y=Basic, color = "Basic"), lwd=1.5) + 
        geom_line(aes(x=x, y=MR, color = "Multiresolution"), lwd=1.5) + 
        theme(legend.title = element_blank()) + 
        theme(legend.position="bottom")
}

