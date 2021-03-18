library(targets)
library(tarchetypes)

options(
  tidyverse.quiet = TRUE,
  clustermq.scheduler = "multicore"
)

# targets:
#
# figures 1 & 2. relative performance of estimators under model1
# figure 3. convergence
# figure 4. relative performance of estimators under model2

tar_option_set(
  packages = c(
    "fastRG",
    "Matrix",
    "RSpectra",
    "cowplot",
    "glue",
    "tidyverse",
    "here",
    "fastadi"
  )
)

##### helper functions ---------------------------------------------------------

# model1


##### end helper functions -----------------------------------------------------

target_estimates <- tar_map(

  unlist = FALSE,

  values = tibble::tibble(
    estimator = rlang::syms(
      c(
        "full_observed_svd" ,
        "adaptive_impute",
        "soft_impute",
        "soft_impute_rank",
        "soft_impute_als",
        "soft_impute_als_rank"
      )
    )
  ),

  tar_target(
    estimate,
    purrr::map(sample, ~estimator(.x, rank = parameters$rank)),
    pattern = map(sample, parameters)
  ),

  tar_target(
    loss,
    purrr::map_dfr(estimate, ~performance(population, .x, params = parameters)),
    pattern = map(population, estimate, parameters)
  )
)

target_combined <- tar_combine(
  combined_losses,
  target_estimates[2],
  command = dplyr::bind_rows(!!!.x, .id = "estimator")
)

num_reps <- 5

list(

  ### reproduce figures 1 & 2

  tar_target(n, 1700),

  tar_target(d, 1000),

  tar_target(sigma, c(0.1, 1, 5, 10, 25, 50)),

  tar_target(p, c(0.05, 0.10, 0.15, 0.20, 0.25)),

  tar_target(rank, c(5, 10, 20, 50)),

  tar_target(
    parameters,
    tibble::tibble(n = n, d = d, sigma = sigma, p = p, rank = rank),
    pattern = cross(n, d, sigma, p, rank)
  ),

  tar_target(
    population,
    model(
      n = parameters$n,
      d = parameters$d,
      sigma = parameters$sigma,
      p = parameters$p,
      rank = parameters$rank
    ),
    pattern = map(parameters),
    iteration = "list"
  ),

  tar_target(
    sample,
    purrr::map(1:num_reps, ~sample_sparse(population)),
    pattern = map(population),
    iteration = "list"
  ),

  target_estimates,

  target_combined,

  tar_target(
    figure1,
    make_figure1(combined_losses)
  ),

  tar_target(
    figure2,
    make_figure1(combined_losses)
  )

  ### reproduce figure 4
)

