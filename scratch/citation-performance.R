# stress test for adaptive-impute-citation performance

library(fastadi)
library(Matrix)
library(profvis)

set.seed(27)

n <- 100000
nnz <- 10000
r <- 10

M <- rsparsematrix(nrow = n, ncol = n, nnz = nnz)
pryr::object_size(M)

s <- RSpectra::svds(M, r)

lra <- citation_adaptive_impute(M, r, max_iter = 5, epsilon = 1e-2)
