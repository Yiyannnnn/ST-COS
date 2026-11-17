seed_text <- Sys.getenv("STCOS_METHOD_SEED", unset = "")
if (nzchar(seed_text)) {
  seed <- suppressWarnings(as.integer(seed_text))
  if (is.na(seed) || seed < 0L) {
    stop("STCOS_METHOD_SEED must be a nonnegative integer.")
  }
  set.seed(seed)
  rm(seed)
}
rm(seed_text)
