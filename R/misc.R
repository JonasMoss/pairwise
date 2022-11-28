#' Calculate asymptotic covariance matrix from singular covariane
#'
#' @param mat Singular covariance matrix.
#' @param i Index of the baseline variable.
#' @param keep_i If `TRUE`, returns a singular covariance matrix with `0` in
#'   the appropriate position.
#' @return Covariance matrix with said baseline variable.

sing_to_cov <- function(mat, i, keep_i = FALSE) {
  k <- nrow(mat)
  mat_i <- mat[i, ] %*% t(rep(1, k))
  temp <- mat - mat_i - t(mat_i) + mat[i, i] * matrix(1, k, k)
  indices <- setdiff(seq(k), i)
  if (keep_i) temp else temp[indices, indices]
}

#' Matrix of resistance distances.
#' @param mat Moore-Penrose inverse of the Laplacian matrix.
#' @return The resistance distance matrix R.
r_matrix <- \(mat) {
  out <- mat * 0
  for (i in seq(nrow(mat))) {
    for (j in seq(nrow(mat))) {
      out[i, j] <- mat[i, i] + mat[j, j] - 2 * mat[i, j]
    }
  }
  out
}

#' Replace Infinite with 0.
#' @param x Vector of inputs.
#' @return Vector with `0`s replacing `Inf`s.
inf_to_zero <- function(x) {
  x[is.infinite(x)] <- 0
  x
}

#' Variance function of binary regression model
#' @param x Vector of inputs.
#' @return Variance at that vector.
varf <- \(x) x * (1 - x)


#' Construct w.
w_new <- function(d_new, beta0, sigma0) {
  mu0 <- c(pnorm(d_new %*% beta0 / sigma0))
  dmu <- c(dnorm(d_new %*% beta0 / sigma0) / sigma0)
  vf <- mu0 * (1 - mu0)
  inf_to_zero(dmu^2 / vf)
}

#' Calculate influence of an observation.
#'
#' Estimates the Kirchhoff coefficient when new observations are added to
#'    the data matrix.
#'
#' @param obj Fitted `pairwise` object.
#' @param source,target The source and target nodes.
#' @param sigma The idiosyncratic standard deviation.
#' @param binary If `TRUE`, the question is binary. Continuous otherwise.
#' @param normalized If `TRUE` returns the old Kirchhoff index divided by the
#'   old Kirchhoff index. If `FALSE`, returns the new Kirchhoff index.
#' @return New normalized or non-normalized Kirchhoff index.
influence <- function(obj, source, target, sigma = 1, binary = FALSE, normalized = TRUE) {
  d_bin <- obj$d_bin
  d_cont <- obj$d_cont
  d_new <- rep(0, ncol(d_cont))
  d_new[c(source, target)] <- c(1, -1) / sigma

  if (binary) {
    w <- c(obj$w, w_new(d_new, obj$beta, obj$sigma))
    d_bin <- rbind(d_bin, d_new)
  } else {
    w <- obj$w
    d_cont <- rbind(d_cont, d_new)
  }

  j <- (crossprod(d_cont) + crossprod(d_bin, d_bin * w)) / obj$sigma^2
  r <- r_matrix(MASS::ginv(j))

  if (normalized) 2 * obj$kirchhoff / sum(r) else 0.5 * sum(r)
}
