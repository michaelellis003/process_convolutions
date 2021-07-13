make_basis <- function(x, knots, radius, form = "uniform") {
    # evaluate the kernel/basis in 1d
    if (!(form %in% c("uniform", "gaussian", "exponential", "polynomial"))) {
        stop('form must be either "uniform" or "gaussian" or "exponential" or "polynomial" ...')
    }
    n       <- length(x)
    n_knots <- length(knots)
    basis   <- matrix(0, n, n_knots)
    for (i in 1:n_knots) {
        if (form == "uniform") {
            basis[, i] <- ifelse(abs(x - knots[i]) < radius, 1, 0)
        }
        if (form == "gaussian") {
            basis[, i] <- dnorm((x - knots[i]) / radius)
        }
        if (form == "polynomial") {
            basis[, i] <- ifelse(abs(x - knots[i]) < radius, (1 - abs(x - knots[i])^3 / radius^3)^3, 0)
        }
        if (form == "exponential") {
            basis[, i] <- dexp(abs((x - knots[i]) / radius))
        }
    }
    return(basis)
}


N <- 1000
x <- seq(0, 1, length.out = N)
n_knots <- 500
knots <- seq(0, 1, length.out = n_knots)
Z <- make_basis(x, knots, radius = 0.2, form = "gaussian")

# plot kernels
matplot(x, Z, type = 'l')

alpha <- rnorm(n_knots)


# white noise
plot(knots, alpha, type = 'h')
abline(h = 0)

# take the convolution
plot(x, Z %*% alpha, type = 'l')





##
## try a gaussian kernel
##

Z <- make_basis(x, knots, radius = 0.2, form = "gaussian")
matplot(x, Z, type = 'l')

# take the convolution
plot(x, Z %*% alpha, type = 'l')

##
## try an exponential kernel
##

Z <- make_basis(x, knots, radius = 0.01, form = "exponential")
matplot(x, Z, type = 'l')

# take the convolution
plot(x, Z %*% alpha, type = 'l')


##
## try a compact suppoted kernel
##

Z <- make_basis(x, knots, radius = 0.1, form = "polynomial")
matplot(x, Z, type = 'l')

# take the convolution
plot(x, Z %*% alpha, type = 'l')
mean(Z == 0)




# Fit to data

y <- rnorm(N) + Z %*% alpha
plot(x, Z %*% alpha, type = 'l')
points(x, y)

n <- 500
s <- sample(1:N, n)

fit <- lm(y[s] ~ Z[s, ] - 1)
alpha_fit <- coefficients(fit)
# plot the fitted line
plot(x[s], y[s])
lines(x, Z %*% alpha_fit, col = "red")


## smoothenss
n_knots <- 1000
knots <- seq(-0.5, 1.5, length.out = n_knots)
Z <- make_basis(x, knots, radius = 0.1, form = "polynomial")
alpha <- rnorm(n_knots)
layout(matrix(1:4, 2, 2))
plot(x, Z %*% alpha, type = 'l', ylab = "Zalpha")
plot(alpha, type = 'l', ylab = "alphas")

# using an autoregressive structure on alpha
alpha[1] <- rnorm(1)
for (i in 2:n_knots) {
    alpha[i] <- rnorm(1) + alpha[i-1]
}
plot(x, Z %*% alpha, type = 'l', ylab = "Zalpha")
plot(alpha, type = 'l', ylab = "alphas")

