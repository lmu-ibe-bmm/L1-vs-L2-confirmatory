# Refresh the agreed study plan and export without fitting models.

update_planning <- function() {
  scripts <- c("run/10_summarize_planning.R", "run/11_export_results.R")
  for (script in scripts) {
    sys.source(script, envir = new.env(parent = globalenv()))
  }
  return(invisible(TRUE))
}
update_planning()
