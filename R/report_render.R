#' Render MethylkeyReport to HTML/PDF
#'
#' Render a complete MethylkeyReport to HTML, PDF, or other formats
#' using Quarto.
#'
#' @param report A \code{MethylkeyReport} object (or character path to report dir).
#' @param output_format Character. Output format: "html", "pdf", "docx", or "all"
#'   (default: "html").
#' @param execute Logical. If TRUE, execute code chunks (default: TRUE).
#'   Set FALSE to skip computation and use cached results.
#' @param quiet Logical. If TRUE, suppress rendering messages (default: FALSE).
#'
#' @return Invisibly returns the report object or path.
#'
#' @details
#' This function:
#' 1. Validates report structure (checks files exist)
#' 2. Generates or updates _quarto.yml if needed
#' 3. Calls quarto::quarto_render()
#' 4. Reports success/failure
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
#' render_report(report, output_format = "html")
#' }
render_report <- function(
    report,
    output_format = "html",
    execute = TRUE,
    quiet = FALSE
) {

  if (is.character(report)) {
    report_dir <- report
  } else if (methods::is(report, "MethylkeyReport")) {
    report_dir <- get_report_subdir(report)
  } else {
    stop("report must be a MethylkeyReport object or character path")
  }

  # Validate structure
  if (!dir.exists(report_dir)) {
    stop(sprintf("Report directory not found: %s", report_dir))
  }

  qmd_files <- list.files(report_dir, pattern = "\\.qmd$")
  if (length(qmd_files) == 0) {
    stop("No .qmd files found in report directory")
  }

  if (!quiet) {
    message(sprintf("\n🚀 Rendering MethylkeyReport: %s", basename(report_dir)))
    message(sprintf("   Files to render: %s", paste(qmd_files, collapse = ", ")))
  }

  # Ensure the report homepage exists for book rendering
  index_path <- file.path(report_dir, "index.qmd")
  if (!file.exists(index_path)) {
    if (methods::is(report, "MethylkeyReport")) {
      report <- create_report_index(report)
    } else {
      stop("Cannot generate index.qmd without report metadata")
    }
  }

  # Generate or update _quarto.yml
  quarto_yml_path <- file.path(report_dir, "_quarto.yml")
  if (!file.exists(quarto_yml_path)) {
    if (methods::is(report, "MethylkeyReport")) {
      generate_quarto_yml(report)
    } else {
      stop("Cannot generate _quarto.yml without report metadata")
    }
  }

  # Call quarto render
  tryCatch({
    if (!requireNamespace("quarto", quietly = TRUE)) {
      stop("quarto package required. Install with: install.packages('quarto')")
    }

    quarto::quarto_render(
      input = report_dir,
      output_format = output_format,
      execute = execute,
      execute_dir = report_dir,
      quiet = quiet
    )

    if (!quiet) {
      message("✓ Report rendered successfully!")
      output_html <- file.path(report_dir, "_output", "index.html")
      if (file.exists(output_html)) {
        message(sprintf("   Output: %s", output_html))
      }
    }

  }, error = function(e) {
    message(sprintf("✗ Rendering failed: %s", e$message))
    stop(e)
  })

  invisible(report)
}

#' Generate _quarto.yml Configuration
#'
#' Create or update the _quarto.yml file for a MethylkeyReport based on
#' generated .qmd files.
#'
#' @param report A \code{MethylkeyReport} object.
#'
#' @return Invisibly returns path to generated _quarto.yml.
#'
#' @details
#' The generated _quarto.yml includes:
#' - Project metadata (title, author, date)
#' - Chapters ordered by .qmd filename (01_, 02_, etc.)
#' - HTML output configuration
#' - Code execution options
#'
#' @keywords internal
generate_quarto_yml <- function(report) {

  if (!methods::is(report, "MethylkeyReport")) {
    stop("report must be a MethylkeyReport object")
  }

  report_dir <- get_report_subdir(report)

  # Get .qmd files in order
  qmd_files <- list.files(report_dir, pattern = "\\.qmd$") |>
    sort()

  if (length(qmd_files) == 0) {
    stop("No .qmd files found in report")
  }

  # Exclude generated pages that should not become book chapters.
  qmd_chapters <- setdiff(qmd_files, c("index.qmd", "model.qmd"))

  # Build chapter list
  chapters <- if (length(qmd_chapters) > 0) {
    paste0("    - ", qmd_chapters, collapse = "\n")
  } else {
    ""
  }

  # Generate YAML content
  yaml_content <- sprintf(
    "project:
  type: book

book:
  title: \"%s\"
  author: \"%s\"
  date: today
  chapters:
    - index.qmd
%s

format:
  html:
    theme: cosmo
    toc: true
    toc-depth: 3
    page-layout: full
    code-fold: true
    code-summary: \"Show code\"
    code-tools: true
    lightbox: true

execute:
  warning: false
  error: false
  freeze: true
",
    report@title,
    report@author,
    chapters
  )

  # Write _quarto.yml
  yaml_path <- file.path(report_dir, "_quarto.yml")
  writeLines(yaml_content, con = yaml_path)

  message(sprintf("✓ Generated _quarto.yml (%d chapters)", length(qmd_chapters)))

  invisible(yaml_path)
}

#' Get Report Output Files
#'
#' List all rendered output files for a report.
#'
#' @param report A \code{MethylkeyReport} object.
#'
#' @return Character vector of output file paths.
#'
#' @export
get_report_outputs <- function(report) {

  if (!methods::is(report, "MethylkeyReport")) {
    stop("report must be a MethylkeyReport object")
  }

  output_dir <- report@output_dir

  if (!dir.exists(output_dir)) {
    return(character(0))
  }

  list.files(output_dir, recursive = TRUE, full.names = TRUE)
}

#' Open Report in Browser
#'
#' Open the rendered HTML report in the default browser.
#'
#' @param report A \code{MethylkeyReport} object.
#'
#' @return Invisibly returns the path to index.html (or NULL if not found).
#'
#' @export
open_report <- function(report) {

  if (!methods::is(report, "MethylkeyReport")) {
    stop("report must be a MethylkeyReport object")
  }

  report_dir <- get_report_subdir(report)
  index_path <- file.path(report_dir, "_output", "index.html")

  if (!file.exists(index_path)) {
    warning(sprintf(
      "Report output not found: %s\nHave you run render_report() yet?",
      index_path
    ))
    return(invisible(NULL))
  }

  message(sprintf("Opening: %s", index_path))

  # Open in browser (platform-dependent)
  if (.Platform$OS.type == "windows") {
    shell.exec(index_path)
  } else if (Sys.info()["sysname"] == "Darwin") {
    system(sprintf("open '%s'", index_path))
  } else {
    system(sprintf("xdg-open '%s'", index_path))
  }

  invisible(index_path)
}

#' Remove Report Project
#'
#' Delete an entire report project directory (use with caution).
#'
#' @param report A \code{MethylkeyReport} object (or character path).
#' @param confirm Logical. If TRUE, require user confirmation (default: TRUE).
#'
#' @return Invisibly returns TRUE if deleted, FALSE otherwise.
#'
#' @export
delete_report <- function(report, confirm = TRUE) {

  if (is.character(report)) {
    report_dir <- report
  } else if (methods::is(report, "MethylkeyReport")) {
    report_dir <- get_report_subdir(report)
  } else {
    stop("report must be a MethylkeyReport object or character path")
  }

  if (!dir.exists(report_dir)) {
    warning(sprintf("Directory not found: %s", report_dir))
    return(invisible(FALSE))
  }

  if (confirm) {
    response <- readline(
      sprintf("Delete report directory? %s (yes/no): ", report_dir)
    )
    if (!tolower(response) %in% c("yes", "y")) {
      message("Cancelled.")
      return(invisible(FALSE))
    }
  }

  unlink(report_dir, recursive = TRUE)
  message(sprintf("✓ Deleted: %s", report_dir))

  invisible(TRUE)
}
