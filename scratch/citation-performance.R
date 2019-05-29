# stress test for adaptive-impute-citation performance

library(fastadi)
library(Matrix)
library(profvis)

set.seed(27)

n <- 1000
nnz <- 1000

M <- rsparsematrix(nrow = n, ncol = n, nnz = nnz)
pryr::object_size(M)

lra <- citation_adaptive_impute(M, 10, epsilon = 1e-4)

profvis({
  lra <- citation_adaptive_impute(M, 10, epsilon = 1e-4)
})

n2 <- 100000
M2 <- matrix(0, nrow = n2, ncol = n2)
pryr::object_size(M2)
