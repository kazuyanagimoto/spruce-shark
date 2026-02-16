tar_model <- tar_plan(
  tar_file_read(
    gg09,
    here_rel("output", "gg09.yaml"),
    yaml::read_yaml(!!.x)
  )
)
