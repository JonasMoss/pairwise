#' Simulate incidence matrix.
#'
#' @param clusters Number of clusters. Defaults to `3`.
#' @param vertices Number of vertices per cluster. Defaults to `8`.
#' @param conn Symmetric matrix with number of vertices between each cluster.
#' @return Simulated data.

simulation <- function(clusters, vertices, conn) {
  d <- matrix(
    0,
    nrow = sum(conn[lower.tri(conn, diag = TRUE)]),
    ncol = clusters * vertices
  )

  replacements <- \(x) rbind(
    cbind(seq(x), seq(x)),
    arrangements::combinations(x, 2)
  )

  indices <- replacements(clusters)
  combinations <- arrangements::combinations(vertices, 2)
  with_replacement <- replacements(vertices)

  from <- 1
  n <- nrow(combinations)
  for (row in seq(nrow(indices))) {
    size <- conn[indices[row, 1], indices[row, 2]]
    if (indices[row, 1] == indices[row, 2]) {
      n <- nrow(combinations)
      i <- combinations[sample(n, size, replace = TRUE), ] +
        matrix(rep((indices[row, ] - 1) * vertices, size), ncol = 2, byrow = TRUE)
    } else {
      n <- nrow(with_replacement)
      i <- with_replacement[sample(n, size, replace = TRUE), ] +
        matrix(rep((indices[row, ] - 1) * vertices, size), ncol = 2, byrow = TRUE)
    }
    i <- matrix(i, ncol = 2)
    elems <- seq(from = from, length.out = size)
    for (j in seq(nrow(i))) {
      d[elems[j], i[j, 1]] <- -1
      d[elems[j], i[j, 2]] <- 1
    }
    from <- from + size
  }

  d * (2 * stats::rbinom(nrow(d), 1, 0.5) - 1)
}
