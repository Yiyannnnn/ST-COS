
library(data.table)
get_process_memory <- function() {
  pid <- 458056
  cmd <- paste0("ps -o rss= -p ", pid)
  mem_kb <- suppressWarnings(as.numeric(system(cmd, intern = TRUE)))
  if (is.na(mem_kb)) {
    cmd <- paste0("ps -p ", pid, " -o rss --no-headers")
    mem_kb <- suppressWarnings(as.numeric(system(cmd, intern = TRUE)))
  }
  if (is.na(mem_kb)) return(0)
  return(mem_kb / 1024)
}

peak_memory <- 0
while (!file.exists("./Reference_based/comp_eff_test/stop_monitor_stcos_50.flag")) {
  current_memory <- get_process_memory()
  if (current_memory > peak_memory) {
    peak_memory <- current_memory
  }
  Sys.sleep(0.5)
}
writeLines(as.character(peak_memory), "./Reference_based/comp_eff_test/peak_memory_stcos_50.txt")

