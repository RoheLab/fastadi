library(RSpectra)

initial_estimate <- function(M, r) {
  p_hat <- mean(M != 0)

  MtM <- crossprod(M)
  MMt <- M %*% t(M)

  sigma_p <- MtM - (1 - p_hat) * diag(MtM)
  sigma_t <- MMt - (1 - p_hat) * diag(MMt)

  svd_p <- svds(sigma_p, r)
  svd_t <- svds(sigma_t, r)

  v_hat <- svd_p$v
  u_hat <- svd_t$u

  n <- nrow(M)
  d <- ncol(M)

  # trace(sigma_p) = sum(diag(sigma_p)) is the sum of eigenvalues of sigma_p

  alpha <- (sum(diag(sigma_p)) - sum(svd_p$d)) / (d - r)
  lambda_hat <- sqrt(svd_p$d - alpha) / p_hat

  svd_M <- svds(M, r)

  # do not understand this yet (asymptotical sign seleciton)
  v_sign <- crossprod(rep(1, d), svd_M$v * v_hat)
  u_sign <- crossprod(rep(1, n), svd_M$u * u_hat)
  s_hat <- c(sign(v_sign * u_sign))

  u_hat %*% diag(s_hat * lambda_hat) %*% t(v_hat)
}

# eigengap heuristic
guestimate_r <- function(M) {
  r <- round(min(dim(M)) / 10)
  s_val <- svds(M, r, nu = 0, nv = 0)$d
  s_ratio <- (s_val[-r] - s_val[-1]) / sval[-r]
  which.max(s_ratio[-1]) + 1
}

adaptive_impute <- function(M, r, eps = 1e-7) {

  d <- ncol(M)
  Z <- initial_estimate(M, r)
  delta <- Inf

  while (delta > eps) {
    M_tilde <- M + Z * (1 - y)  # sanity check

    svd_M <- svds(M_tilde, r)

    u_hat <- svd_M$u
    v_hat <- svd_M$v

    # get all the eigenvalues squared somehow through this?
    alpha <- (sum(colSums(M_tilde^2)) - sum(svd_M$d^2)) / (d - r)

    # will explode if you get a negative here
    lambda_hat <- sqrt(svd_M$d^2 - alpha)

    Z_new <- u_hat %*% diag(lambda_hat) %*% t(v_hat)

    delta <- sum((Z_new - Z)^2) / sum(Z^2)
    Z <- Z_new

    print(glue::glue("Delta: {round(delta, 8)}, alpha: {round(alpha, 3)}"))
  }

  Z
}

n <- 500
d <- 100
r <- 5

A <- matrix(runif(n * r, -5, 5), n, r)
B <- matrix(runif(d * r, -5, 5), d, r)
M0 <- A %*% t(B)

err <- matrix(rnorm(n * d), n, d)
Mf <- M0 + err

p <- 0.3
y <- matrix(rbinom(n * d, 1, p), n, d)
dat <- Mf * y

init <- initial_estimate(dat, r)
filled <- adaptive_impute(dat, r, eps = 1e-2)
filled2 <- adaptive_impute(dat, r, eps = 1e-9)
filled3 <- adaptive_impute(dat, r, eps = 1e-11)

sum((init - dat)^2)
sum((filled - dat)^2)
sum((filled2 - dat)^2)
sum((filled3 - dat)^2)


filled2 <- adaptive_impute(dat, r)
sum((filled2 - dat)^2) / prod(dim(dat))


