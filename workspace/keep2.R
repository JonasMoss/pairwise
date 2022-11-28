set.seed(313)

clusters <- 3
vertices <- 5
cross <- 2
conn <- matrix(c(
  25, 2, 2,
  2, 25, 5,
  2, 5, 25
), nrow = 3, byrow = TRUE)

sigma <- c(
  rep(1, sum(diag(conn))),
  rep(2, (sum(conn) - sum(diag(conn))) / 2)
)

d <- simulate(clusters, vertices, conn)
beta <- seq(ncol(d)) - 1
tau <- 1
cont <- rnorm(nrow(d), d %*% beta, tau * sigma)
bin <- cont >= 0
cont[rbinom(nrow(d), 1, 0.3) == 1] <- NA

data <- data.frame(
  cont = cont,
  bin = bin,
  sigma = sigma,
  d
)
