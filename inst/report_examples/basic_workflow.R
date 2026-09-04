#' Example: Progressive Methylation Analysis Report
#'
#' This script demonstrates how to use methylkey's Phase 6 report builder
#' to construct a complete methylation analysis report step-by-step.
#'
#' The approach is progressive: you can add steps incrementally, modify
#' parameters, and rebuild specific steps without recalculating earlier work.

library(methylkey)

# ============================================================================
# 1. INITIALIZE A NEW REPORT
# ============================================================================

# Create a new report project with a unique ID
# This allows multiple reports in the same folder
report <- create_methylkey_report(
  project_dir = "~/methylation_analysis",
  report_id = "exp_001_raw",
  idat_path = "/path/to/idat/files",
  sample_sheet_path = "/path/to/sample_sheet.csv",
  title = "Methylation Analysis - Experiment 001",
  author = "Your Name"
)

# View report status
print(report)
report_status(report)

# ============================================================================
# 2. ADD ANALYSIS STEPS PROGRESSIVELY
# ============================================================================

# Step 1: Add data loading chapter
report <- report |>
  add_report_data_loading(
    prep = "TQCDPB",
    ncore = 4
  )

# Step 2: Add quality control chapter
report <- report |>
  add_report_qc(
    groups = c("group", "age", "sentrix_id", "gender")
  )

# Step 3: Add exploratory analysis (PCA, MDS)
report <- report |> add_report_group_control(
  groups = c("group", "age", "sentrix_id", "gender"),
  columns = c("samples", "group", "age", "sentrix_id", "gender"),
  sva = NULL,
  win = TRUE,
  xy = FALSE
)

# ============================================================================
# 3. ADD MULTIPLE DIFFERENTIAL ANALYSIS MODELS
# ============================================================================

# Model 1: Simple group comparison
report <- report |>
  add_report_model_analysis(
    model = "~group",
    intercept = "WT",
    method = "ls",
    dmrtools = c("dmrcate", "ipdmr"),
    genome = "hg38"
  )

# Model 2: Group + Age adjustment
report <- report |>
  add_report_model_analysis(
    model = "~group+age",
    intercept = "WT",
    dmrtools = c("dmrcate", "ipdmr", "combp")
  )

# Model 3: Group only in young mice
# (This would require custom filtering - example for advanced usage)
report <- report |>
  add_report_subgroup(
    subgroup_var = "age",
    subgroup_values = c("young"),
    model = "~group",
    intercept = "WT"
  )

# ============================================================================
# 4. RENDER THE COMPLETE REPORT
# ============================================================================

# Check status before rendering
report_status(report)

# Render everything to HTML
render_report(report, output_format = "html")

# Open in browser
open_report(report)

# ============================================================================
# 5. POST-ANALYSIS CLEANUP AND MANIPULATION
# ============================================================================

# If you want to rebuild a specific step (e.g., after changing parameters)
report <- rebuild_report_step(report, step = "qc", force = TRUE)

# Clear Quarto cache to save space after completing analysis
clear_report_cache(report)

# Get list of output files
output_files <- get_report_outputs(report)

# ============================================================================
# ADVANCED: MULTIPLE REPORTS IN SAME FOLDER
# ============================================================================

# You can create separate reports for different analyses in the same directory
# using different report_ids

report_sva <- create_methylkey_report(
  project_dir = "~/methylation_analysis",
  report_id = "exp_001_sva_corrected",  # Different ID
  idat_path = "/path/to/idat/files",
  sample_sheet_path = "/path/to/sample_sheet.csv",
  title = "Methylation Analysis - SVA Corrected"
) |>
  add_report_data_loading(...) |>
  add_report_qc(...) |>
  add_report_group_control(correction = "sva", ...)  # Different correction

render_report(report_sva)

# Both reports coexist without conflict:
# ~/methylation_analysis/
#   ├── exp_001_raw/
#   │   ├── 01_data_loading.qmd
#   │   ├── ...
#   │   └── _output/
#   └── exp_001_sva_corrected/
#       ├── 01_data_loading.qmd
#       ├── ...
#       └── _output/

# ============================================================================
# WORKFLOW BENEFITS
# ============================================================================

# ✓ Progressive: Add steps incrementally without redoing previous work
# ✓ Modular: Each step is independent
# ✓ Reproducible: Script is complete documentation
# ✓ Flexible: Easy to add/modify/rebuild individual steps
# ✓ Scalable: Support for multiple models, subgroups, etc.
# ✓ Multi-report: Avoid file conflicts when testing different approaches
# ✓ Cached: Quarto handles incremental computation automatically
