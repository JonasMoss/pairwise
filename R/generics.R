#' Plot connection strength graph of a fitted pairwise model.
#'
#' @param x Fitted pairwise model.
#' @param induced If `TRUE`, plots induced connections strengths in addition
#'    to the main connection strengths.
#' @param col Passed to `col`.
#' @param reference Reference node.
#' @param labels Labels of nodes. Only used if `reference` is `NULL`.
#' @param ... Ignored.
#' @export

plot.pairwise <- function(x, reference = NULL, induced = TRUE, col = "black", labels = seq(obj$p), ...) {
  obj = x
  if (!is.null(reference)) {
    y <- stats::coef(obj, reference = reference)
    x <- seq(length(y))
    cis <- stats::confint(obj, reference = reference)
    Hmisc::errbar(
      x = x, y = y, yplus = cis[, 2], yminus = cis[, 1],
      ylab = "Value", xlab = "Question index", type = "b"
    )
    graphics::grid()
    Hmisc::errbar(
      x = x, y = y, yplus = cis[, 2], yminus = cis[, 1],
      add = TRUE
    )
  } else {
    data <- obj$data
    p <- obj$p
    n <- obj$n

    eig <- eigen(obj$j_inv)
    singulars <- eig$values
    singulars[length(singulars)] <- 0
    vecs <- eig$vectors
    res <- umap::umap(t(diag(singulars) %*% t(vecs)))$layout
    plot(res, type = "n", xlab = "x", ylab = "y",
      main = "Model-implied graph"
    )
    r <- r_matrix(obj$j_inv)
    graphics::text(res, labels = labels, col = col)

    combs <- arrangements::combinations(p, 2, replace = FALSE)

    d <- \(left, right) 1 / sqrt(r[left, right])

    if (induced) {
      for (i in seq(n)) {
        graphics::lines(res[combs[i, ], ],
          lwd = d(combs[i, 1], combs[i, 2]), col = "grey",
          lty = 1
        )
      }
    }

    for (i in seq(n)) {
      indices <- which(data[i, -(1:3)] != 0)
      graphics::lines(res[indices, ], lwd = d(indices[1], indices[2]), col = "black")
    }
  }
}

#' @export
confint.pairwise <- function(object, parm, level = 0.95, reference = 1, ...) {
  beta <- stats::coef(object, reference)
  ses <- diag(sing_to_cov(object$j_inv, reference, keep_i = TRUE))
  modifier <- c(ses * stats::qnorm(1 / 2 + level / 2))
  cbind(lower = beta - modifier,upper = beta + modifier)[parm, ]
}

#' @export
coef.pairwise <- function(object, reference = NULL, ...) {
  if (is.null(reference)) object$beta else object$beta - object$beta[reference]
}

#' @export
predict.pairwise <- function(object, sources = NULL, targets = NULL, ...) {
  if (is.null(sources) | is.null(targets)) {
    sources <- which(t(object$data[, -c(1:3)]) == 1, arr.ind = TRUE)[, 1]
    targets <- which(t(object$data[, -c(1:3)]) == -1, arr.ind = TRUE)[, 1]
  }
  stats::coef(object)[sources] - stats::coef(object)[targets]
}
