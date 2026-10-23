# Export public results without any job registries or per-job cache files.

public_files <- function(path) {
  if (!dir.exists(path)) return(character())
  entries <- list.files(path, full.names = TRUE)
  excluded <- grepl("registry|pilot_registry_lock|worker_bootstrap", basename(entries))
  entries <- entries[!excluded]
  files <- entries[!dir.exists(entries)]
  files <- files[grepl("[.](csv|rds|png|pdf|txt|tex|md|sha256)$", files)]
  directories <- entries[dir.exists(entries)]
  return(c(files, unlist(lapply(directories, public_files), use.names = FALSE)))
}
if (!file.exists("outputs/planning_summary/planning_status.txt")) {
  stop("Run the planning summary after both new pilots have completed.")
}

