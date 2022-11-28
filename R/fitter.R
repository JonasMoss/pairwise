#' Fit a pairwise comparison model
#'
#' @param data Data frame to estimate.
#' @param epsilon Congvergence tolerance. Iterations converge when the relative
#'    difference is less than epsilon.
#' @param maxit Maximum number of iterations. Defaults to `25`, which coincides
#'    with that of `glm`.
#' @param beta1 Optional starting value.
#' @export
pairwise <- function(data, maxit = 25, epsilon = 1e-8, start = NULL) {
  stopifnot(is_induced_connected(data))
  bins <- is.na(data$cont)
  d_ind <- 4:ncol(data)
  weights_bin <- 1 / data$sigma[bins]
  weights_cont <- 1 / data$sigma[!bins]

  d_bin <- diag(weights_bin) %*% as.matrix(data[bins, d_ind])
  y_bin <- data[bins, 2]
  d_cont <- diag(weights_cont) %*% as.matrix(data[!bins, d_ind])
  y_cont <- data[!bins, 1] * weights_cont

  n <- nrow(data)
  p <- ncol(data) - 3
  n_cont <- nrow(d_cont)
  div <- max(n_cont - p, 1)

  j <- \(w) crossprod(d_cont) + crossprod(d_bin, d_bin * w)

  beta1 <- if (is.null(start)) rep(0, ncol(data) - 3) else start
  sigma1 <- 1 # Starting value.

  for (i in seq(maxit)) {
    beta0 <- beta1
    sigma0 <- sigma1

    mu <- c(pnorm(d_bin %*% beta0 / sigma0))
    dmu <- c(dnorm(d_bin %*% beta0 / sigma0) / sigma0)
    vf <- mu * (1 - mu)
    d <- inf_to_zero(dmu / vf)
    w <- inf_to_zero(dmu^2 / vf)
    u_bin <- 1 / sigma0 * (t(d_bin) * d) %*% (y_bin - mu)
    u_cont <- 1 / sigma0^2 * t(d_cont) %*% (y_cont - d_cont %*% beta0)
    j_inv <- sigma0^2 * MASS::ginv(j(w))

    beta1 <- c(beta0 + j_inv %*% (u_bin + u_cont))
    sigma1 <- sqrt(sum((d_cont %*% beta1 - y_cont)^2) / div)

    if (sum(abs(beta1 - beta0) / abs(beta0)) < epsilon) break
  }

  if (i == maxit) warning("Maximum number of iterations reached.")

  attr(beta1, "j_inv") <- j_inv
  attr(beta1, "options") <- list(maxit = maxit, epsilon = epsilon, iter = i)
  if (any(w == 0)) attr(beta1, "fitted_0") <- TRUE
  attr(beta1, "w") <- w
  attr(beta1, "sigma") <- sigma1
  attr(beta1, "p") <- p
  attr(beta1, "n") <- n
  attr(beta1, "data") <- data
  class(beta) <- "pairwise"
  beta1
}
