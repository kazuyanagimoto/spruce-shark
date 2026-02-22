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
  tar_figure,
  tar_data,
  tar_model,
  tar_file(
    manuscript_src,
    list.files(
      here_rel("manuscript"),
      pattern = "\\.(qmd|tex|bib|lua)$",
      recursive = TRUE,
      full.names = TRUE
    )
  ),
  tar_file(
    manuscript,
    {
      manuscript_src
      quarto::quarto_render(
        here_rel("manuscript", "spruce_shark_yanagimoto.qmd")
      )
      here_rel("manuscript", "spruce_shark_yanagimoto.pdf")
    }
  )
)
