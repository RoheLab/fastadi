% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/adaptive-impute.R
\name{adaptive_impute}
\alias{adaptive_impute}
\alias{adaptive_impute.sparseMatrix}
\alias{adaptive_impute.LRMF}
\title{AdaptiveImpute}
\usage{
adaptive_impute(
  X,
  rank,
  ...,
  initialization = c("svd", "adaptive-initialize", "approximate"),
  max_iter = 200L,
  check_interval = 1L,
  epsilon = 1e-07,
  additional = NULL
)

\method{adaptive_impute}{sparseMatrix}(
  X,
  rank,
  ...,
  initialization = c("svd", "adaptive-initialize", "approximate"),
  additional = NULL
)

\method{adaptive_impute}{LRMF}(
  X,
  rank,
  ...,
  epsilon = 1e-07,
  max_iter = 200L,
  check_interval = 1L
)
}
\arguments{
\item{X}{A sparse matrix of \code{\link[Matrix:sparseMatrix]{Matrix::sparseMatrix()}} class.}

\item{rank}{Desired rank (integer) to use in the low rank approximation.
Must be at least \code{2L} and at most the rank of \code{X}. Note that the rank
of \code{X} is typically unobserved and computations may be unstable or
even fail when \code{rank} is near or exceeds this threshold.}

\item{...}{Unused additional arguments.}

\item{initialization}{How to initialize the low rank approximation.
Options are:
\itemize{
\item \code{"svd"} (default). In the initialization step, this treats
unobserved values as zeroes.
\item \code{"adaptive-initialize"}. In the initialization step, this treats
unobserved values as actually unobserved. However, the current
\code{AdaptiveInitialize} implementation relies on dense matrix
computations that are only suitable for relatively small matrices.
\item \code{"approximate"}. An approximate variant of \code{AdaptiveInitialize}
that is less computationally expensive. See \code{adaptive_initialize}
for details.
}

Note that initialization matters as \code{AdaptiveImpute} optimizes
a non-convex objective. The current theory shows that initializing
with \code{AdaptiveInitialize} leads to a consistent estimator, but it
isn't know if this is the case for SVD initialization. Empirically
we have found that SVD initialization works well nonetheless.}

\item{max_iter}{Maximum number of iterations to perform (integer). Defaults
to \code{200L}. In practice 10 or so iterations will get you a decent
approximation to use in exploratory analysis, and and 50-100 will get
you most of the way to convergence. Must be at least \code{1L}.}

\item{check_interval}{Integer specifying how often to perform convergence
checks. Defaults to \code{1L}. In practice, check for convergence requires
a norm calculation that is expensive for large matrices and decreasing
the frequency of convergence checks will reduce computation time. Can
also be set to \code{NULL}, which case \code{max_iter} iterations of the algorithm
will occur with no possibility of stopping due to small relative change
in the imputed matrix. In this case \code{delta} will be reported as \code{Inf}.}

\item{epsilon}{Convergence criteria, measured in terms of relative change
in Frobenius norm of the full imputed matrix. Defaults to \code{1e-7}.}

\item{additional}{Ignored except when \code{alpha_method = "approximate"}
in which case it controls the precise of the approximation to \code{alpha}.
The approximate computation of \code{alpha} will always understand \code{alpha},
but the approximation will be better for larger values of \code{additional}.
We recommend making \code{additional} as large as computationally tolerable.}
}
\value{
A low rank matrix factorization represented by an
\code{\link[=adaptive_imputation]{adaptive_imputation()}} object.
}
\description{
An implementation of the \code{AdaptiveImpute} algorithm for matrix completion
for sparse matrices.
}
\examples{

### SVD initialization (default) --------------------------------------------

mf <- adaptive_impute(ml100k, rank = 3L, max_iter = 20L)

# build a rank-5 approximation only for
# observed elements of ml100k

preds <- predict(mf, ml100k)

# estimate the in-sample reconstruction mse

R <- resid(mf, ml100k)
norm(R, type = "F") / nnzero(ml100k)

### Exact AdaptiveInitialize initialization ---------------------------------

mf2 <- adaptive_impute(
  ml100k,
  rank = 3L,
  max_iter = 20L,
  initialization = "adaptive-initialize"
)

R2 <- resid(mf2, ml100k)
norm(R2, type = "F") / nnzero(ml100k)

### Approximate AdaptiveInitialize initialization ---------------------------

mf3 <- adaptive_impute(
  ml100k,
  rank = 3L,
  max_iter = 20L,
  initialization = "approximate",
  additional = 25
)

R3 <- resid(mf3, ml100k)
norm(R3, type = "F") / nnzero(ml100k)

}