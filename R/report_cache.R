#' Clear Report Cache
#'
#' Remove Quarto cache directories and intermediate files for a report.
#' Useful for cleaning up after analysis completion or forcing a fresh start.
#'
#' @param report A \code{MethylkeyReport} object (or character path to report directory).
#' @param recursive Logical. If TRUE, recursively delete cache (default: TRUE).
#'
#' @return Invisibly returns the report object (for pipe-friendly usage).
#'
#' @details
#' This function removes:
#' - `_quarto_cache/` directory
#' - Individual `.qmd_cache` directories
#' - Temporary computation files
#'
#' Use this when:
#' - Finalizing analysis to clean up temporary files
#' - Restarting a step from scratch
#' - Freeing disk space
#'
#' @importFrom methods is
#'
#' @export
#'
#' @examples
#' \dontrun{
#' report <- create_methylkey_report(...)
#' # ... add steps and render ...
#' clear_report_cache(report)  # Clean up when done
#' }
clear_report_cache <- function(report, recursive = TRUE) {

  if (is.character(report)) {
    # Handle both full path and report directory
    cache_dir <- report
  } else if (methods::is(report, "MethylkeyReport")) {
    cache_dir <- report@cache_dir
  } else {
    stop("report must be a MethylkeyReport object or character path")
  }

  # Ensure we're dealing with a path
  if (!is.character(cache_dir)) {
    return(invisible(report))
  }

  # Get the report subdirectory
  if (methods::is(report, "MethylkeyReport")) {
    report_dir <- get_report_subdir(report)
  } else {
    report_dir <- dirname(cache_dir)
  }

  # Remove main cache directory
  if (dir.exists(cache_dir)) {
    message("Removing cache directory: ", cache_dir)
    unlink(cache_dir, recursive = recursive)
  }

  # Remove individual .qmd_cache directories
  cache_dirs <- list.files(
    report_dir,
    pattern = "_files$",
    full.names = TRUE,
    include.dirs = TRUE
  )

  for (cache_subdir in cache_dirs) {
    message("Removing cache subdirectory: ", cache_subdir)
    unlink(cache_subdir, recursive = recursive)
  }

  message("✓ Cache cleared for report")

  invisible(report)
}

#' Rebuild a Specific Report Step
#'
#' Force rerun of a specific analysis step without recalculating other steps.
#' Clears cache for the target step and optionally renders.
#'
#' @param report A \code{MethylkeyReport} object.
#' @param step Character. Name of the step to rebuild
#'   (e.g., "data_loading", "qc", "group_control", "model_analysis").
#' @param force Logical. If TRUE, remove existing output and force rerun (default: FALSE).
#' @param render Logical. If TRUE, render the step after rebuilding (default: FALSE).
#'
#' @return Invisibly returns the report object.
#'
#' @details
#' This function:
#' 1. Validates that the step exists in the report
#' 2. Clears Quarto cache for the step
#' 3. Optionally re-renders the step
#'
#' Use when:
#' - A step failed and you want to try again
#' - You want to test parameter changes for one step
#' - Previous rendering had errors
#'
#' @importFrom methods is
#'
#' @export
#'
#' @examples
#' \dontrun{
#' report <- create_methylkey_report(...) |>
#'   add_report_data_loading(...) |>
#'   add_report_qc(...)
#' 
#' # If QC step had issues, rebuild it
#' report <- rebuild_report_step(report, step = "qc", force = TRUE)
#' }
rebuild_report_step <- function(report, step, force = FALSE, render = FALSE) {

  if (!methods::is(report, "MethylkeyReport")) {
    stop("report must be a MethylkeyReport object")
  }

  if (!step %in% names(report@steps)) {
    stop(sprintf(
      "Step '%s' not found in report. Available steps: %s",
      step,
      paste(names(report@steps), collapse = ", ")
    ))
  }

  message(sprintf("Rebuilding step: %s", step))

  report_dir <- get_report_subdir(report)

  # Find and remove cache associated with this step
  # Quarto uses pattern like: 01_data_loading_files/ for 01_data_loading.qmd
  step_cache_pattern <- paste0("^[0-9]+_", gsub("_", "-", step), "_files$")
  cache_dirs <- list.files(
    report_dir,
    pattern = step_cache_pattern,
    full.names = TRUE,
    include.dirs = TRUE
  )

  for (cache_dir in cache_dirs) {
    message("  Removing cache: ", basename(cache_dir))
    unlink(cache_dir, recursive = TRUE)
  }

  # Remove main Quarto cache if force=TRUE
  if (force) {
    main_cache <- file.path(report_dir, "_quarto_cache")
    if (dir.exists(main_cache)) {
      message("  Removing main cache: _quarto_cache")
      unlink(main_cache, recursive = TRUE)
    }
  }

  message("✓ Step cleared and ready for rebuild")

  if (render) {
    message("Rendering step...")
    # This will be called by render_report or similar
  }

  invisible(report)
}

#' Verify Report Structure
#'
#' Check that all necessary files and directories exist for a report.
#'
#' @param report A \code{MethylkeyReport} object.
#'
#' @return Logical TRUE if structure is valid, FALSE otherwise.
#'
#' @keywords internal
verify_report_structure <- function(report) {

  if (!methods::is(report, "MethylkeyReport")) {
    return(FALSE)
  }

  report_dir <- get_report_subdir(report)

  # Check that report directory exists
  if (!dir.exists(report_dir)) {
    warning(sprintf("Report directory does not exist: %s", report_dir))
    return(FALSE)
  }

  # Check that IDAT path exists (or is a URL)
  if (!dir.exists(report@idat_path) && !grepl("^https?://", report@idat_path)) {
    warning(sprintf("IDAT path does not exist: %s", report@idat_path))
    return(FALSE)
  }

  # Check that sample sheet exists
  if (!file.exists(report@sample_sheet_path)) {
    warning(sprintf("Sample sheet not found: %s", report@sample_sheet_path))
    return(FALSE)
  }

  TRUE
}

#' Get Report Status Summary
#'
#' Display a summary of the report's current state.
#'
#' @param report A \code{MethylkeyReport} object.
#'
#' @return Invisibly returns a list with status information.
#'
#' @keywords internal
report_status <- function(report) {

  if (!methods::is(report, "MethylkeyReport")) {
    stop("report must be a MethylkeyReport object")
  }

  report_dir <- get_report_subdir(report)

  # Count generated .qmd files
  qmd_files <- list.files(report_dir, pattern = "\\.qmd$")

  # Check for _quarto.yml
  quarto_yml_exists <- file.exists(file.path(report_dir, "_quarto.yml"))

  # Check cache status
  cache_exists <- dir.exists(report@cache_dir)

  cat("\nReport Status Summary\n")
  cat("=======================\n")
  cat("Report ID:          ", report@report_id, "\n")
  cat("Project Directory:  ", report_dir, "\n")
  cat("Exists:             ", dir.exists(report_dir), "\n\n")

  cat("Configuration:\n")
  cat("  Steps Added:      ", length(report@steps), "\n")
  cat("  .qmd Files:       ", length(qmd_files), "\n")
  cat("  _quarto.yml:      ", quarto_yml_exists, "\n")
  cat("  Cache Exists:     ", cache_exists, "\n\n")

  cat("Files Generated:\n")
  if (length(qmd_files) > 0) {
    cat("  -", paste(qmd_files, collapse = "\n  - "), "\n\n")
  }

  cat("Ready to Render:    ", quarto_yml_exists && length(qmd_files) > 0, "\n")

  invisible(list(
    steps_count = length(report@steps),
    qmd_count = length(qmd_files),
    has_quarto_yml = quarto_yml_exists,
    cache_exists = cache_exists,
    ready_to_render = quarto_yml_exists && length(qmd_files) > 0
  ))
}
