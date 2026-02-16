here_rel <- function(...) {
  fs::path_rel(here::here(...))
}

download_file <- function(url, destfile, ...) {
  if (!file.exists(destfile)) {
    download.file(url, destfile, ...)
  }
  return(destfile)
}

write_lines <- function(text, file, ...) {
  writeLines(text = text, con = file, ...)
  return(file)
}
