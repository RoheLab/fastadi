##### the following is from Adi_Rfunction (original) and Adi_manual (original)

#####################
## Adaptive-Impute ##
#####################

#' Adaptive-Impute
#'
#' An iterative matrix completion algorithm based on the thresholded SVD.
#' It differentially and adaptively penalizes the singular valuesin each
#' iteration. Although Adaptive-Impute employs multiple thresholding
#' parameters updated every iteration, there is no tuning problem since it
#' automatically finds specific values of the thresholding parameters which
#' are theoretically-justified and data-dependent.
#'
#' For more details, please see “Intelligent Initialization and Adaptive
#' Thresholding for Iterative Matrix Completion; Some Statistical and
#' Algorithmic Theory for Adaptive-Impute”.
#'
#' @param M.p input matrix where the unobserved matrix are marked as zeros.
#'  Can be  in sparse matrix format (inherit from class "sparseMatrix"
#'  as in package Matrix)
#' @param r rank of the resulting matrix. If no value is given, it
#'  automatically chooses one with the biggest eigengap based on a scree plot
#' @param sparse `TRUE` if only a few entries are observed and `FALSE` otherwise
#' @param sign.choice we recommend to use “asympt”. It uses asymptotically
#'  consistent signs when combining the estimates of singular values
#'  and singular vectors to obtain the estimates of the low-rank matrix.
#'  “lm” uses linear regression and “greedy” uses greedy search.
#' @param min.value a constant. lower bound for the entries
#' @param max.value a constant. upper bound for the entries
#' @param tol tolerance for convergence
#' @param itmax maximum number of iterations
#'
#' @return TODO
#' @export
#'
#' @author Juhee Cho, Donggyu Kim, Karl Rohe
#'
#' @examples
#'
#' n <- 500
#' d <- 100
#' r <- 5
#'
#' A <- matrix(runif(n * r, -5, 5), n, r)
#' B <- matrix(runif(d * r, -5, 5), d, r)
#' M0 <- A %*% t(B)
#'
#' err <- matrix(rnorm(n * d), n, d)
#' Mf <- M0 + err
#'
#' p <- 0.1
#' y <- matrix(rbinom(n * d, 1, p), n, d)
#' dat <- Mf * y
#'
#' out <- adaptImpute(M.p = dat, r = r)
adaptImpute <- function(M.p, r = NULL, sparse = c(TRUE, FALSE),
                        sign.choice = c("asympt", "lm", "greedy"),
                        min.value = NULL, max.value = NULL,
                        tol = 1e-07, itmax = 200) {
  n <- nrow(M.p)
  d <- ncol(M.p)
  ind <- !(M.p == 0)
  if (is.null(r)) {
    r.temp <- round(min(n, d) / 10)
    sval <- svds(M.p, 94, nu = 0, nv = 0)$d
    sval.ratio <- (sval[1:(r.temp - 1)] - sval[2:r.temp]) / sval[1:(r.temp - 1)]
    r <- which.max(sval.ratio[-1]) + 1
  }
  p.hat <- mean(ind)
  temp.M <- Initial(M.p, r, p.hat, n, d, ind, sign.choice)
  if (sparse == TRUE | p.hat < 0.5) {
    M.p <- Matrix(M.p, sparse = TRUE)
    ifelse(is.null(max.value), max.value <- max(M.p@x), temp.M[temp.M > max.value] <- max.value)
    ifelse(is.null(min.value), min.value <- min(M.p@x), temp.M[temp.M < min.value] <- min.value)
  } else {
    ifelse(is.null(max.value), max.value <- max(M.p[ind]), temp.M[temp.M > max.value] <- max.value)
    ifelse(is.null(min.value), min.value <- min(M.p[ind]), temp.M[temp.M < min.value] <- min.value)
  }

  itr <- 0
  error <- Inf
  while (error > tol) {
    M <- temp.M * (1 - ind) + M.p
    temp <- SVD.F(M, r, M.p, d)
    temp.M1 <- temp
    error <- sum((temp.M1 - temp.M)^2) / sum(temp.M^2)
    temp.M <- temp.M1
    max.ind <- temp.M > max.value
    temp.M <- temp.M * (1 - max.ind) + max.value * max.ind
    min.ind <- temp.M < min.value
    temp.M <- temp.M * (1 - min.ind) + min.value * min.ind
    itr <- itr + 1
    if (itr %% 10 == 0) cat(".")
    if (itr > itmax) {
      break
    }
  }
  return(temp.M)
}

#####################################
## subfunctions for adaptiveImpute ##
#####################################

sign.lm <- function(M.p, lambda.hat, u.hat, v.hat, r, ind) {
  nz <- which(ind)
  X <- matrix(0, sum(ind), r)
  for (i in 1:r) {
    E <- lambda.hat[i] * u.hat[, i] %*% t(v.hat[, i])
    X[, i] <- E[nz]
  }
  y <- M.p[nz]
  XX.inv <- solve(crossprod(X))
  XTY <- t(X) %*% y
  coef <- XX.inv %*% XTY
  lambda.hat.R <- diag(c(lambda.hat * coef))
  out <- u.hat %*% lambda.hat.R %*% t(v.hat)
  return(out)
}

sign.asympt <- function(M.p, n, d, lambda.hat, u.hat, v.hat, r, ind) {
  svd.Mp <- svds(M.p, r)
  Vp.Vhat <- crossprod(rep(1, d), (svd.Mp$v * v.hat))
  Up.Uhat <- crossprod(rep(1, n), (svd.Mp$u * u.hat))
  coef <- sign(Vp.Vhat * Up.Uhat)
  lambda.hat.R <- diag(c(lambda.hat * coef))
  out <- u.hat %*% (lambda.hat.R %*% t(v.hat))
  return(out)
}

sign.greedy <- function(M.p, lambda.hat, u.hat, v.hat, r, ind) {
  nr <- 2^r
  coef.set <- matrix(0, nr, r)
  for (i in 1:r) { # i=1
    motiv <- c(rep(1, nr / (2^i)), rep(-1, nr / (2^i)))
    coef.set[, i] <- rep(motiv, 2^(i - 1))
  }
  MSEs <- c()
  for (j in 1:nr) { # j=1
    coef <- coef.set[j, ]
    lambda.hat.R <- diag(c(lambda.hat * coef))
    Mhat <- u.hat %*% (lambda.hat.R %*% t(v.hat))
    diff <- (M.p - Mhat)[ind]
    MSEs[j] <- mean(diff^2)
  }
  pick <- which.min(MSEs)
  coef <- coef.set[pick, ]
  lambda.hat.R <- diag(c(lambda.hat * coef))
  out <- u.hat %*% (lambda.hat.R %*% t(v.hat))
  return(out)
}

Initial <- function(M.p, r, p.hat, n, d, ind, sign.choice) {
  SIG <- (t(M.p) %*% M.p) / p.hat^2
  SIG.diag <- (p.hat - 1) * diag(SIG)
  diag(SIG) <- diag(SIG) + SIG.diag
  SIG2 <- (M.p %*% t(M.p)) / p.hat^2
  SIG2.diag <- (p.hat - 1) * diag(SIG2)
  diag(SIG2) <- diag(SIG2) + SIG2.diag

  obj <- svds(SIG, r)
  obj2 <- svds(SIG2, r)

  tau <- (sum(diag(SIG)) - sum(obj$d)) / (d - r)
  lambda.hat <- sqrt(obj$d - tau)
  v.hat <- obj$v
  u.hat <- obj2$v

  if (sign.choice == "asympt") {
    M.hat <- sign.asympt(M.p, n, d, lambda.hat, u.hat, v.hat, r, ind)
  } else {
    if (sign.choice == "greedy") {
      M.hat <- sign.greedy(M.p, lambda.hat, u.hat, v.hat, r, ind)
    } else {
      M.hat <- sign.lm(M.p, lambda.hat, u.hat, v.hat, r, ind)
    }
  }
  return(M.hat)
}

SVD.F <- function(M, r, M.p, d) {
  obj <- svds(M, r)
  M2 <- M^2
  MTM <- colSums(M2)
  tau <- (sum(MTM) - sum(obj$d^2)) / (d - r)
  lambda.hat <- sqrt(obj$d^2 - tau)
  v.hat <- obj$v
  u.hat <- obj$u
  M.hat <- u.hat %*% (diag(lambda.hat) %*% t(v.hat))
  return(M.hat)
}
