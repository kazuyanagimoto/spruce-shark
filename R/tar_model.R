tar_model <- tar_plan(
  tar_map(
    values = lst(name = c("model", "estimation", "benchmark")),
    tar_file_read(
      jl,
      here_rel("Julia", paste0(name, ".jl")),
      readLines(!!.x)
    )
  ),
  tar_file(
    output,
    run_jl(
      jl_file_estimation,
      c(
        here_rel("output", "ct_estimated.yaml"),
        here_rel("output", "moment.csv"),
        here_rel("output", "simulation.csv")
      ),
      jl_file_model
    ),
  ),
  tar_file_read(ct, output[1], yaml::read_yaml(!!.x)),
  tar_file_read(moment, output[2], readr::read_csv(!!.x)),
  tar_file_read(simulation, output[3], readr::read_csv(!!.x)),
  tar_file_read(
    benchmark,
    run_jl(
      jl_file_benchmark,
      here_rel("output", "benchmark.csv"),
      jl_file_model
    ),
    readr::read_csv(!!.x)
  )
)

run_jl <- function(jl_file, output_files, ...) {
  cmd <- sprintf(
    'julia --project="%s" "%s"',
    here::here(),
    jl_file
  )
  status <- system(cmd)
  if (status != 0) {
    stop("Julia script failed: ", jl_file)
  }
  return(output_files)
}
