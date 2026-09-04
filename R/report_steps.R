#' Template Rendering Utilities
#'
#' Internal functions for loading and rendering .qmd templates with parameter substitution.
#'
#' @keywords internal

#' Load Template File
#'
#' Read a template file from inst/templates/
#'
#' @param template_name Character. Name of template file (with .template extension).
#'
#' @return Character string with template content.
#'
#' @keywords internal
load_template <- function(template_name) {
  template_path <- system.file(
    "templates",
    template_name,
    package = "methylkey"
  )

  if (!file.exists(template_path)) {
    template_path <- file.path("inst", "templates", template_name)
  }

  if (!file.exists(template_path)) {
    stop(sprintf("Template not found: %s", template_name))
  }

  readLines(template_path, warn = FALSE) |>
    paste(collapse = "\n")
}

#' Render Template with Parameters
#'
#' Replace placeholders in template with actual values using glue.
#'
#' @param template Character. Template content with `{{PLACEHOLDER}}` markers.
#' @param params List. Named list of values to interpolate.
#'
#' @return Character string with rendered content.
#'
#' @importFrom glue glue_data
#'
#' @keywords internal
render_template <- function(template, params) {
  params <- lapply(params, function(value) {
    if (is.null(value)) "NULL" else value
  })

  glue::glue_data(params, template, .open = "{{", .close = "}}") |>
    as.character()
}

#' Create Report Home Page
#'
#' Generate the report's `index.qmd` home page with analysis metadata and
#' session information.
#'
#' @param report A \code{MethylkeyReport} object.
#'
#' @return The updated \code{MethylkeyReport} object.
#'
#' @keywords internal
create_report_index <- function(report) {
  if (!methods::is(report, "MethylkeyReport")) {
    stop("report must be a MethylkeyReport object")
  }

  template <- load_template("index.qmd.template")
  params <- list(
    TITLE = report@title,
    AUTHOR = report@author,
    DATE = report@date,
    REPORT_ID = report@report_id,
    PROJECT_DIR = report@project_dir,
    IDAT_PATH = report@idat_path,
    SAMPLE_SHEET_PATH = report@sample_sheet_path,
    CREATED_AT = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  )

  content <- render_template(template, params)
  write_qmd_file(report, "index.qmd", content)

  report@created_files <- c(report@created_files, "index.qmd")
  report@metadata$last_modified <- Sys.time()

  invisible(report)
}

#' Write QMD File
#'
#' Write rendered .qmd content to report directory.
#'
#' @param report A \code{MethylkeyReport} object.
#' @param filename Character. Name of output .qmd file.
#' @param content Character. Rendered template content.
#'
#' @return Character path to written file (invisibly).
#'
#' @keywords internal
write_qmd_file <- function(report, filename, content) {
  report_dir <- get_report_subdir(report)
  output_path <- file.path(report_dir, filename)

  writeLines(content, con = output_path)
  message(sprintf("  Generated: %s", filename))

  invisible(output_path)
}

#' Copy Model Template
#'
#' Add the reusable model template to a report when model analysis is used.
#'
#' @param report A \code{MethylkeyReport} object.
#'
#' @return Invisibly returns the report object.
#'
#' @keywords internal
copy_model_template <- function(report) {
  report_dir <- get_report_subdir(report)
  model_path <- file.path(report_dir, "model.qmd")

  if (!file.exists(model_path)) {
    write_qmd_file(report, "model.qmd", load_template("model.qmd"))
    report@created_files <- unique(c(report@created_files, "model.qmd"))
  }

  invisible(report)
}

#' Add Data Loading Step
#'
#' Generate the data loading chapter (.qmd) for the report.
#'
#' @param report A \code{MethylkeyReport} object.
#' @param prep Character. sesame preprocessing method (default: "TQCDPB").
#' @param ncore Integer. Number of cores for parallel processing (default: 4).
#' @param clock_models List. Named list of clock model paths (optional).
#'
#' @return The updated \code{MethylkeyReport} object (invisibly, for piping).
#'
#' @export
#'
#' @examples
#' \dontrun{
#' report <- create_methylkey_report(...) |>
#'   add_report_data_loading(prep = "TQCDPB", ncore = 4)
#' }
add_report_data_loading <- function(
    report,
    prep = "TQCDPB",
    ncore = 4,
    clock_models = NULL
) {

  if (!methods::is(report, "MethylkeyReport")) {
    stop("report must be a MethylkeyReport object")
  }

  message("Adding step: data_loading")

  # Load and render template
  template <- load_template("01_data_loading.qmd.template")

  # Prepare parameters for interpolation
  params <- list(
    IDAT_PATH = report@idat_path,
    SAMPLE_SHEET_PATH = report@sample_sheet_path,
    PREP = prep,
    NCORE = ncore,
    DEBUG = TRUE
  )

  if (!is.null(clock_models)) {
    # Format clock_models for template
    clock_str <- paste(
      names(clock_models),
      "=",
      paste0("'", unname(clock_models), "'"),
      collapse = ",\n    "
    )
    params$CLOCK_MODELS <- clock_str
  } else {
    params$CLOCK_MODELS <- ""
  }

  content <- render_template(template, params)
  write_qmd_file(report, "01_data_loading.qmd", content)

  # Update report object
  report@steps$data_loading <- list(
    template = "01_data_loading.qmd.template",
    params = list(prep = prep, ncore = ncore, clock_models = clock_models)
  )
  report@created_files <- c(report@created_files, "01_data_loading.qmd")
  report@metadata$last_modified <- Sys.time()
  generate_quarto_yml(report)

  invisible(report)
}

#' Add Quality Control Step
#'
#' Generate the QC chapter (.qmd) for the report.
#'
#' @param report A \code{MethylkeyReport} object.
#' @param groups Character vector. Variables to visualize in QC plots
#'   (default: c("sentrix_id", "group", "age", "gender")).
#'
#' @return The updated \code{MethylkeyReport} object (invisibly, for piping).
#'
#' @export
add_report_qc <- function(
  report,
  groups = c("sentrix_id", "group", "age", "gender")
) {

  if (!methods::is(report, "MethylkeyReport")) {
    stop("report must be a MethylkeyReport object")
  }

  if ("data_loading" %notin% names(report@steps)) {
    warning("data_loading step not added yet. QC depends on it.")
  }

  message("Adding step: qc")

  template <- load_template("02_quality_control.qmd.template")

  # Format groups parameter
  groups_r <- paste0("c('", paste(groups, collapse = "', '"), "')")

  params <- list(
    GROUPS = groups_r
  )

  content <- render_template(template, params)
  write_qmd_file(report, "02_quality_control.qmd", content)

  report@steps$qc <- list(
    template = "02_quality_control.qmd.template",
    params = list(groups = groups)
  )
  report@created_files <- c(report@created_files, "02_quality_control.qmd")
  report@metadata$last_modified <- Sys.time()
  generate_quarto_yml(report)

  invisible(report)
}

#' Add Group Control Step
#'
#' Generate the group control / exploratory analysis chapter (.qmd).
#'
#' @param report A \code{MethylkeyReport} object.
#' @param groups Character vector. Grouping variables to analyze.
#' @param sva Character. SVA model path (optional).
#' @param win Logical. Apply winsorization (default: TRUE).
#' @param xy Logical. Exclude sex chromosomes (default: FALSE).
#'
#' @return The updated \code{MethylkeyReport} object (invisibly, for piping).
#'
#' @export
add_report_group_control <- function(
    report,
    groups = c("sentrix_id", "group", "age", "gender"),
    columns = c("sample_id", "sentrix_id", "group", "age", "gender"),
    sva = NULL,
    win = TRUE,
    xy = FALSE
) {

  if (!methods::is(report, "MethylkeyReport")) {
    stop("report must be a MethylkeyReport object")
  }

  message("Adding step: group_control")

  template <- load_template("03_group_control.qmd.template")

  groups_r <- paste0("c('", paste(groups, collapse = "', '"), "')")
  columns_r <- paste0("c('", paste(columns, collapse = "', '"), "')")

  params <- list(
    GROUPS = groups_r,
    COLUMNS = columns_r,
    SVA = sva,
    WIN = win,
    XY = xy
  )

  content <- render_template(template, params)
  write_qmd_file(report, "03_group_control.qmd", content)

  report@steps$group_control <- list(
    template = "03_group_control.qmd.template",
    params = list(
      groups = groups,
      sva = sva,
      win = win,
      xy = xy
    )
  )
  report@created_files <- c(report@created_files, "03_group_control.qmd")
  report@metadata$last_modified <- Sys.time()
  generate_quarto_yml(report)

  invisible(report)
}

#' Add Model Analysis Step
#'
#' Generate a differential analysis chapter for a specific model.
#' Can be called multiple times to add multiple models.
#'
#' @param report A \code{MethylkeyReport} object.
#' @param model Character. Model formula (e.g., "~group" or "~group+age").
#' @param intercept Character. Reference level for the first grouping variable.
#' @param method Character. Fitting method: "ls" or "robust" (default: "ls").
#' @param dmrtools Character vector. DMR tools to run
#'   (default: c("dmrcate", "ipdmr")).
#' @param genome Character. Genome assembly (default: "hg38").
#' @param model_id Character. Optional unique identifier for this model
#'   (default: auto-generated from model string).
#'
#' @return The updated \code{MethylkeyReport} object (invisibly, for piping).
#'
#' @export
#'
#' @examples
#' \dontrun{
#' report <- create_methylkey_report(...) |>
#'   add_report_model_analysis(
#'     model = "~group",
#'     intercept = "WT",
#'     dmrtools = c("dmrcate", "ipdmr")
#'   ) |>
#'   add_report_model_analysis(
#'     model = "~group+age",
#'     intercept = "WT"
#'   )
#' }
add_report_model_analysis <- function(
    report,
    model,
    intercept,
    method = "ls",
    dmrtools = c("dmrcate", "ipdmr"),
    genome = "hg38",
    model_id = NULL
) {

  if (!methods::is(report, "MethylkeyReport")) {
    stop("report must be a MethylkeyReport object")
  }

  # Auto-generate model_id if not provided
  if (is.null(model_id) || model_id == "") {
    model_id <- make.names(tolower(gsub("[^a-zA-Z0-9]", "_", model)))
    model_id <- gsub("^_+|_+$", "", model_id)
  }

  message(sprintf("Adding step: model_analysis (%s)", model_id))

  template <- load_template("04_model_analysis.qmd.template")

  dmrtools_r <- paste0("c('", paste(dmrtools, collapse = "', '"), "')")

  params <- list(
    MODEL = model,
    INTERCEPT = intercept,
    METHOD = method,
    DMRTOOLS = dmrtools_r,
    GENOME = genome,
    MODEL_ID = model_id,
    TXDB = case_when(
      genome == "hg38" ~ "TxDb.Hsapiens.UCSC.hg38.knownGene",
      genome == "hg19" ~ "TxDb.Hsapiens.UCSC.hg19.knownGene",
      genome == "mm10" ~ "TxDb.Mmusculus.UCSC.mm10.knownGene",
      TRUE ~ ""
    )
  )

  content <- render_template(template, params)

  # Filename based on model_id
  filename <- sprintf("04_model_%s.qmd", model_id)
  write_qmd_file(report, filename, content)

  # Store with unique key (handle multiple models)
  step_key <- sprintf("model_analysis_%s", model_id)
  report@steps[[step_key]] <- list(
    template = "04_model_analysis.qmd.template",
    params = list(
      model = model,
      intercept = intercept,
      method = method,
      dmrtools = dmrtools,
      genome = genome
    ),
    model_id = model_id
  )
  report@created_files <- c(report@created_files, filename)
  report@metadata$last_modified <- Sys.time()
  report <- copy_model_template(report)
  generate_quarto_yml(report)

  invisible(report)
}

#' Add Subgroup Analysis Step
#'
#' Generate a subgroup-specific analysis chapter.
#'
#' @param report A \code{MethylkeyReport} object.
#' @param subgroup_var Character. Variable name to split on (e.g., "age_group").
#' @param subgroup_values Character vector. Values to include
#'   (if NULL, uses all unique values).
#' @param model Character. Model formula for subgroup analysis.
#' @param intercept Character. Reference level.
#'
#' @return The updated \code{MethylkeyReport} object (invisibly, for piping).
#'
#' @export
add_report_subgroup <- function(
    report,
    subgroup_var,
    subgroup_values = NULL,
    model = "~group",
    intercept = "WT"
) {

  if (!methods::is(report, "MethylkeyReport")) {
    stop("report must be a MethylkeyReport object")
  }

  message(sprintf("Adding step: subgroup_analysis (%s)", subgroup_var))

  template <- load_template("05_subgroup_analysis.qmd.template")

  if (!is.null(subgroup_values)) {
    subgroup_values_r <- paste0("c('", paste(subgroup_values, collapse = "', '"), "')")
  } else {
    subgroup_values_r <- "NULL"
  }

  params <- list(
    SUBGROUP_VAR = subgroup_var,
    SUBGROUP_VALUES = subgroup_values_r,
    MODEL = model,
    INTERCEPT = intercept
  )

  content <- render_template(template, params)

  filename <- sprintf("05_subgroup_%s.qmd", tolower(gsub(" ", "_", subgroup_var)))
  write_qmd_file(report, filename, content)

  step_key <- sprintf("subgroup_%s", subgroup_var)
  report@steps[[step_key]] <- list(
    template = "05_subgroup_analysis.qmd.template",
    params = list(
      subgroup_var = subgroup_var,
      subgroup_values = subgroup_values,
      model = model,
      intercept = intercept
    )
  )
  report@created_files <- c(report@created_files, filename)
  report@metadata$last_modified <- Sys.time()
  generate_quarto_yml(report)

  invisible(report)
}

#' Helper: %notin% operator
#'
#' Negation of %in%
#'
#' @keywords internal
`%notin%` <- function(x, y) !(x %in% y)
