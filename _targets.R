library(targets)
library(tarchetypes)
suppressPackageStartupMessages(library(dplyr))

options(
  dplyr.summarise.inform = FALSE,
  readr.show_col_types = FALSE
)

tar_config_set(
  store = here::here("_targets"),
  script = here::here("_targets.R")
)

tar_source()
tar_plan(
  tar_data,
  tar_model
)
