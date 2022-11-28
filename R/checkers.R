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
directed <- function(d) {
  sources <- which(d == 1, arr.ind = TRUE)
  sources <- sources[order(sources[, 1]), 2]
  targets <- which(d == -1, arr.ind = TRUE)
  targets <- targets[order(targets[, 1]), 2]
  igraph::graph_from_edgelist(cbind(sources, targets))
}
