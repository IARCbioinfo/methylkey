#' MethylkeyReport S4 Class
#'
#' An S4 class representing a methylation analysis report project.
#' Tracks project structure, steps, cache configuration, and metadata.
#'
#' @slot project_dir Character. Root directory for the report project.
#' @slot report_id Character. Unique identifier for this report (allows multiple
#'   reports in same folder).
#' @slot title Character. Report title.
#' @slot author Character. Report author(s).
#' @slot date Character. Report date (default: today).
#' @slot steps List. Named list of added steps with their parameters.
#' @slot idat_path Character. Path to IDAT files.
#' @slot sample_sheet_path Character. Path to sample sheet.
#' @slot cache_dir Character. Directory for Quarto cache files.
#' @slot output_dir Character. Directory for rendered output.
#' @slot created_files Character vector. Paths to generated .qmd files.
#' @slot metadata List. Additional metadata key-value pairs.
#'
#' @name MethylkeyReport
#' @aliases MethylkeyReport-class
#' @docType class
#' @exportClass MethylkeyReport
setClass("MethylkeyReport",
  slots = list(
    project_dir = "character",
    report_id = "character",
    title = "character",
    author = "character",
    date = "character",
    steps = "list",
    idat_path = "character",
    sample_sheet_path = "character",
    cache_dir = "character",
    output_dir = "character",
    created_files = "character",
    metadata = "list"
  ),
  prototype = list(
    project_dir = character(0),
    report_id = character(0),
    title = character(0),
    author = character(0),
    date = character(0),
    steps = list(),
    idat_path = character(0),
    sample_sheet_path = character(0),
    cache_dir = character(0),
    output_dir = character(0),
    created_files = character(0),
    metadata = list()
  )
)

#' Create a New Methylation Analysis Report Project
#'
#' Factory function to initialize a new MethylkeyReport object and set up
#' the project directory structure.
#'
#' @param project_dir Character. Root directory for the report.
#' @param report_id Character. Unique identifier for this report.
#'   Default: random UUID-like string.
#' @param idat_path Character. Path to directory containing IDAT files.
#' @param sample_sheet_path Character. Path to sample sheet CSV/TSV.
#' @param title Character. Report title (default: "Methylation Analysis").
#' @param author Character. Author name(s) (default: "").
#' @param date Character. Date for report (default: today's date).
#'
#' @return An object of class \code{MethylkeyReport}.
#'
#' @importFrom methods new
#'
#' @export
#'
#' @examples
#' \dontrun{
#' report <- create_methylkey_report(
#'   project_dir = "~/my_analysis",
#'   report_id = "exp_001",
#'   idat_path = "/data/idats/",
#'   sample_sheet_path = "samples.csv"
#' )
#' }
create_methylkey_report <- function(
    project_dir,
    report_id = NULL,
    idat_path,
    sample_sheet_path,
    title = "Methylation Analysis",
    author = "",
    date = NULL
) {

  # Validate inputs
  assertthat::assert_that(
    is.character(project_dir) && length(project_dir) == 1,
    msg = "project_dir must be a single character string"
  )
  assertthat::assert_that(
    is.character(idat_path) && length(idat_path) == 1,
    msg = "idat_path must be a single character string"
  )
  assertthat::assert_that(
    is.character(sample_sheet_path) && length(sample_sheet_path) == 1,
    msg = "sample_sheet_path must be a single character string"
  )

  if (is.null(report_id) || report_id == "") {
    report_id <- paste0("report_", format(Sys.time(), "%Y%m%d_%H%M%S"))
  }

  if (is.null(date)) {
    date <- format(Sys.Date(), "%Y-%m-%d")
  }

  # Normalize paths
  project_dir <- normalizePath(project_dir, mustWork = FALSE)
  idat_path <- normalizePath(idat_path, mustWork = FALSE)
  sample_sheet_path <- normalizePath(sample_sheet_path, mustWork = FALSE)

  # Create project directory structure
  if (!dir.exists(project_dir)) {
    dir.create(project_dir, recursive = TRUE)
  }

  # Subdirectories based on report_id
  report_subdir <- file.path(project_dir, report_id)
  if (!dir.exists(report_subdir)) {
    dir.create(report_subdir, recursive = TRUE)
  }

  cache_dir <- file.path(report_subdir, "_quarto_cache")
  output_dir <- file.path(report_subdir, "_output")

  # Create output directories
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Initialize MethylkeyReport object
  report <- methods::new(
    "MethylkeyReport",
    project_dir = project_dir,
    report_id = report_id,
    title = title,
    author = author,
    date = date,
    steps = list(),
    idat_path = idat_path,
    sample_sheet_path = sample_sheet_path,
    cache_dir = cache_dir,
    output_dir = output_dir,
    created_files = character(0),
    metadata = list(
      created_at = Sys.time(),
      last_modified = Sys.time()
    )
  )

  # Create the report homepage so book rendering has a valid home page.
  report <- create_report_index(report)

  # Generate _quarto.yml immediately so the report is ready to render.
  generate_quarto_yml(report)

  message(sprintf(
    "✓ MethylkeyReport '%s' initialized in %s\n",
    report_id,
    report_subdir
  ))

  report
}

#' Get Report Subdirectory
#'
#' Returns the full path to the report's subdirectory.
#'
#' @param report A \code{MethylkeyReport} object.
#'
#' @return Character string with full path.
#'
#' @keywords internal
get_report_subdir <- function(report) {
  file.path(report@project_dir, report@report_id)
}

#' Show Method for MethylkeyReport
#'
#' Display summary information about a MethylkeyReport object.
#'
#' @param object A \code{MethylkeyReport} object.
#'
#' @importFrom methods setMethod
#' @export
setMethod("show", "MethylkeyReport", function(object) {
  cat("MethylkeyReport Object\n")
  cat("=======================\n\n")
  cat("Report ID:        ", object@report_id, "\n")
  cat("Project Dir:      ", object@project_dir, "\n")
  cat("Title:            ", object@title, "\n")
  cat("Author:           ", object@author, "\n")
  cat("Date:             ", object@date, "\n\n")

  cat("Data Configuration:\n")
  cat("  IDAT Path:        ", object@idat_path, "\n")
  cat("  Sample Sheet:     ", object@sample_sheet_path, "\n\n")

  cat("Steps Added:      ", length(object@steps), "\n")
  if (length(object@steps) > 0) {
    cat("  -", paste(names(object@steps), collapse = "\n  - "), "\n\n")
  }

  cat("Files Generated:  ", length(object@created_files), "\n")
  cat("Cache Directory:  ", object@cache_dir, "\n")
  cat("Output Directory: ", object@output_dir, "\n")
})

