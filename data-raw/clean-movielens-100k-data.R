library(readr)
library(dplyr)
library(here)
library(Matrix)

ml100k_df <- read_table2(
  here("data-raw/ml-100k/u.data"),
  col_names = c("user_id", "item_id", "rating", "timestamp")
) %>%
  select(-timestamp) %>%
  mutate_all(as.integer)

ml100k <- sparseMatrix(
  i = ml100k_df$user_id,
  j = ml100k_df$item_id,
  x = ml100k_df$rating
)

rownames(ml100k) <- paste0("user", 1:nrow(ml100k))
colnames(ml100k) <- paste0("item", 1:ncol(ml100k))

usethis::use_data(ml100k, overwrite = TRUE)
