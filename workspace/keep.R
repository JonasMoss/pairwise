set.seed(313)

## First simulate some nodes.
clusters <- 3
vertices <- 5
cross <- 2
conn <- matrix(c(
  25, 2, 2,
  2, 25, 5,
  2, 5, 25
), nrow = 3, byrow = TRUE)

conn = conn

sigma <- c(
  rep(1, sum(diag(conn))),
  rep(2, (sum(conn) - sum(diag(conn)))/2))

d <- simulation(clusters, vertices, conn)
beta <- runif(ncol(d), 0, 10)
tau <- 1
cont <- rnorm(nrow(d), d %*% beta, tau * sigma)
bin <- cont >= 0
cont[rbinom(nrow(d), 1, 0.3) == 1] <- NA

data <- data.frame(
  cont = cont,
  bin = bin,
  sigma = sigma,
  d)







data <- make_frame(dat.research)
data[1, 1] = NA
data$cont = data$cont - 0.11
obj <- pairwise(data)








plot(induced(data))


obj = pairwise(data)
influence(obj, 1, 2, binary = FALSE)
plot.pairwise(
  obj,
  col = c(rep("black", 5,), rep("red", 5), rep("blue", 5)),
  induced = FALSE)

i = 1
plot(obj, reference = i, col = c(rep("black", 5,), rep("red", 5), rep("blue", 5)))
lines(seq(length(beta)), beta - beta[i], lty = 2, lwd = 2)


source <- 5
target <- 6
sigma_new = 1
d_new = rep(0, ncol(data) - 3)
d_new[c(source, target)] = c(1, -1)
y_new <- rnorm(1, d_new %*% beta, tau * sigma_new)

data2 <- rbind(
  data,
  c(y_new, y_new > 0, sigma_new, d_new))


obj2 = pairwise(data2)
plot(obj2, reference = i)


res_cont <- outer(1:15, 1:15, Vectorize(function(i,j) influence(obj2, i, j, binary = FALSE)))
res_bin <- outer(1:15, 1:15, Vectorize(function(i,j) influence(obj2, i, j, binary = TRUE)))

rcont <- round(-res_cont / sum(r_matrix(attr(obj2, "j_inv"))) + 1, 3)
rbin <- round(-res_bin / sum(r_matrix(attr(obj2, "j_inv"))) + 1, 3)


indices <- which(rcont * 0 == 0, arr.ind = TRUE)
cbind(indices[order(rcont, decreasing = TRUE), ], sort(rcont, decreasing = TRUE))

n <- 15
mat <- lower.tri(matrix(0, nrow = n, ncol = n))
indices <- which(mat, arr.ind = TRUE)
base <- sum(r_matrix(attr(obj2, "j_inv")))
for (row in seq(nrow(indices))) {
  i = indices[row, 1]
  j = indices[row, 2]
  mat[i, j] <- round(1-influence(obj, i, j, binary = FALSE)/base)
}

tri <- lower.tri(matrix(0, nrow = n, ncol = n))

cbind(indices[order(mat[tri], decreasing = TRUE), ],
      sort(mat[tri], decreasing = TRUE))


w <- MASS::ginv(attr(obj, "w"))
1 - influence(obj, 5, 6, binary = FALSE)/sum(r_matrix(attr(obj, "j_inv")))
res_cont <- outer(1:15, 1:15, Vectorize(function(i,j) influence(obj, i, j, binary = FALSE)))
res_bin <- outer(1:15, 1:15, Vectorize(function(i,j) influence(obj, i, j, binary = TRUE)))

round(sum(r_matrix(attr(obj, "j_inv"))) / res_cont - 1, 2)
round(sum(r_matrix(attr(obj, "j_inv"))) / res_bin - 1, 2)


influence(obj, 1, 2, binary = FALSE)

j = 1
i = setdiff(1:15, 1)
plot(i, sapply(i, \(i) influence(obj, j, i, binary = TRUE)
               /influence(obj, j, i, binary = FALSE)))

cont <- outer(1:15, 1:15, Vectorize(\(i,j) influence(obj, j, i, binary = TRUE)))
lattice::levelplot(cont)

binary <- outer(1:15, 1:15, Vectorize(\(i,j) influence(obj, j, i, binary = FALSE)))
lattice::levelplot(binary)

res <- outer(1:15, 1:15, Vectorize(\(i,j) influence(obj, j, i, binary = FALSE)))
lattice::levelplot(binary / cont)

which.max(cont)
