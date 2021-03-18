model <- function(
  A, B,
  n, d, r,
  sigma,
  p
) {
  object <- list(
    A = matrix()
    n = n,
    d = d,
    r = r,
    sigma = sigma,
    p = p
  )

  class(object) <- "noisy_low_rank_model"
  object
}

?sample
