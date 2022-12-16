#' Makes a suitable data drame out of pairwise agreement data.
#'
#' @param data Pairwise agreement data.
#' @param keep_names If `TRUE`, keeps the names of the questions. This can
#'   make the data frame unwieldly. Defaults to `FALSE`.
#' @return An appropriate data frame.
make_frame <- function(data, keep_names = FALSE) {
  levels <- levels(as.factor(c(data$source, data$target)))
  source <- as.numeric(factor(data$source, levels = levels))
  target <- as.numeric(factor(data$target, levels = levels))

  k <- length(levels)
  n <- nrow(data)

  d <- matrix(data = 0, nrow = n, ncol = k)
  for (i in seq(n)) {
    d[i, source[i]] <- -1
    d[i, target[i]] <- 1
  }

  selection <- seq(k)
  y <- log(as.numeric(data$distance))

  data_frame <- data.frame(cont = y, bin = TRUE, sigma = 1, d[, selection])

  data_frame
}


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
#' @param d_new New data.
#' @param beta0 Vector of regression coefficients.
#' @param sigma0 The standard deviation.
w_new <- function(d_new, beta0, sigma0) {
  mu0 <- c(stats::pnorm(d_new %*% beta0 / sigma0))
  dmu <- c(stats::dnorm(d_new %*% beta0 / sigma0) / sigma0)
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

#' Check if the induced graph is  connected.
#'
#' @param data A data frame of pairwise comparison data.
#' @param mode Passed to `igraph::is.connected`. Defaults to strong.
#' @return `TRUE` if the induced graph is strongly connected, `FALSE` if not.
is_induced_connected <- function(data, mode = "strong") {
  g <- induced(data)
  igraph::is.connected(g, mode)
}

#' Transform data frame to induced graph.
#'
#' @param data data frame of pairwise comparison data.
#' @return The induced graph.
induced <- function(data) {
  bins <- is.na(data$cont)
  d_ind <- 4:ncol(data)
  graph_bin <- directed(data[bins, d_ind] * (2 * data$bin[bins] - 1))
  graph_cont_1 <- directed(data[!bins, d_ind])
  graph_cont_2 <- directed(-data[!bins, d_ind])
  igraph::union(graph_bin, graph_cont_1, graph_cont_2)
}

#' Construct directed graph from d.
#' @param d Data source.
#' @return A directed graph.
directed <- function(d) {
  sources <- which(d == 1, arr.ind = TRUE)
  sources <- sources[order(sources[, 1]), 2]
  targets <- which(d == -1, arr.ind = TRUE)
  targets <- targets[order(targets[, 1]), 2]
  igraph::graph_from_edgelist(cbind(sources, targets))
}
