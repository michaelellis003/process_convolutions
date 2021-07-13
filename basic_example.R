basic_model <- function() {
    library(ggplot2)
    library(gridExtra)
    library(latex2exp)
    source("make_kernel.R")
    source("mcmc_basic.R")
    
    ## set ggplot theme
    theme_set(theme_minimal())
    
    ## Simulate data ---------------------------------------------------------------
    N <- 50
    x <- seq(0, 1, length.out = N)
    n_knots <- 10
    knots <- seq(0, 1, length.out = n_knots)
    
    sigma2_epsilon <- 0.1^2
    sigma2_alpha <- 1
    alpha <- rnorm(n_knots, 0, sqrt(sigma2_alpha))
    
    Z <- make_kernel(x, knots, stddev = 0.25)
    
    y <- Z %*% alpha + rnorm(N, 0, sqrt(sigma2_epsilon))
    
    ## basic model -----------------------------------------------------------------
    n_mcmc <- 5000
    burnin <- 2500
    basic_fit <- mcmc_basic(y, x, 
                            knots, stddev=0.25, 
                            n_mcmc = n_mcmc, burnin = burnin, n_message = 500)
    
    plot_mcmc <- data.frame(
        Samples = 1:(n_mcmc-burnin),
        sigma2_epsilon = basic_fit$sigma2_epsilon,
        sigma2_alpha = basic_fit$sigma2_alpha
    )
    
    p1 <- ggplot() + 
        geom_line(data = plot_mcmc, aes(x=Samples, y=sigma2_epsilon)) + 
        geom_hline(yintercept = sigma2_epsilon, color = "red") +
        ylab(TeX(r'($\sigma^2_{\epsilon}$)'))
    
    p2 <- ggplot() + 
        geom_line(data = plot_mcmc, aes(x=Samples, y=sigma2_alpha)) + 
        geom_hline(yintercept = sigma2_alpha, color = "red") +
        ylab(TeX(r'($\sigma^2_{\alpha}$)'))
    
    plot_data <- data.frame(
        x = x,
        z = Z %*% alpha,
        y = y,
        Basic = Z %*% apply(basic_fit$alpha, 2, mean)
    )
    
    p3 <- ggplot(data = plot_data) + 
        geom_point(aes(x=x, y=y)) + 
        geom_line(aes(x=x, y=z, color = "True")) + 
        geom_line(aes(x=x, y=Basic, color = "Basic")) +
        theme(legend.title = element_blank()) + 
        theme(legend.position="bottom")
    
    grid.arrange(p1, p2, p3, nrow = 3)
}