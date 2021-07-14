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