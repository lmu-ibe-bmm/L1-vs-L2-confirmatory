# Select this file with source(file.choose()) after extracting the project ZIP.
# Locate the project and refresh its completed planning results.

start_project <- function() {
  frames <- sys.frames()
  scripts <- unlist(lapply(frames, function(frame) {
    if (is.null(frame$ofile)) character() else as.character(frame$ofile)
  }), use.names = FALSE)
  scripts <- scripts[basename(scripts) == "start.R"]
  if (!length(scripts)) {
    stop("Open ridge-lasso-comparison and run source('start.R').")
  }
  root <- dirname(normalizePath(tail(scripts, 1L), mustWork = TRUE))
  required <- c("run_config.R", "R/pilot_core.R", "update_planning.R")
  if (!all(file.exists(file.path(root, required)))) {
    stop("The project folder is incomplete. Extract the whole ZIP first.")
  }
  setwd(root)
  message("Project folder: ", getwd())
  sys.source("update_planning.R", envir = new.env(parent = globalenv()))
  return(invisible(root))
}
start_project()
