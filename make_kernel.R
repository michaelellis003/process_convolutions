make_kernel <- function(x, knots, stddev) {
    
    n <- length(x)
    n_knots <- length(knots)
    d <- fields::rdist(x, knots)
    kern <- matrix(0, n, n_knots)
    kern <- dnorm(d, 0, stddev)
    
    return(kern)
}