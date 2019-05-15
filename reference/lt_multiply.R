n <- 7
d <- 5
k <- 3

# want to calculate: y = (L * U V') x  // where * is hadamard prod
#
# L: n x d
# U V': n x d
#
# k is the rank of decomp
#
# U: n x k
# V: d x k
#
# x: d-dimensional vector we would like to multiply by

M <- matrix(rep(0, n * d), n, d)
L <- M
L[lower.tri(L)] <- 1
L

U <- matrix(1:(n * k), n, k)
V <- matrix(rev(1:(d * k)), d, k)

x <- 1:d

#

# sanity check

stopifnot(dim(L) == c(n, d))
stopifnot(dim(U %*% t(V)) == c(n, d))

# brute force approach -- H for Hadamard product (elementwise)

H <- (L * U %*% t(V))
y <- H %*% x
y

# clever (hopefully approach)

f <- function() {

  y <- rep(0, n)

  for (i in 1:n) {
    for (j in 1:d) {
      y[i] <- y[i] + L[i, j] * sum(U[i, ] * V[j, ]) * x[j]
    }
  }

  y
}

stopifnot(f() == y)

# clever 2

f2 <- function() {

  y <- rep(0, n)

  for (i in 1:n) {
    for (j in 1:d) {
      t <- 0
      for (l in 1:k) {
        t <- t + U[i, l] * V[j, l]
      }
      y[i] <- y[i] + L[i, j] * t * x[j]
    }
  }

  y
}
f2()
stopifnot(f2() == y)


# clever 3

f3 <- function() {

  y <- rep(0, n)

  for (i in 1:n) {
    t <- 0
    for (l in 1:k) {
      tt <- 0
      for (j in 1:d) {
        tt <- tt + L[i, j] *  V[j, l] * x[j]
      }
      t <- t + U[i, l] * tt
    }
    y[i] <- t
  }

  y
}

f3()
stopifnot(f3() == y)

## i think the question is how much of the computation can be vectorized


# sanity check: dense L too large for memory



