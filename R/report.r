#' Generate Interactive Methylation Report
#'
#' Creates an interactive Quarto report with Shiny widgets for dynamic exploration
#' of differential methylation analysis results. Users can switch between models,
#' statistical methods (ls/robust), data types (raw/SVA-corrected), and visualization tools.
#'
#' @param se SummarizedExperiment object with raw beta values
#' @param mrs MethylResultSet object with differential methylation results
#' @param mvals_raw Matrix of raw M-values
#' @param mvals_sva Matrix of SVA-corrected M-values
#' @param genome Genome build ("hg19", "hg38", "mm10", etc.)
#' @param dmrtools Character vector of DMR detection tools to include
#' @param output_file Path where to save the .qmd file (default: "methylkey_report.qmd")
#' @param title Report title
#' @param author Report author name
#'
#' @return Invisibly returns the path to the generated .qmd file
#'
#' @details
#' The generated report includes:
#' - Interactive model selector (for multiple contrasts)
#' - Method selector (ls vs robust regression)
#' - Data type selector (raw vs SVA-corrected M-values)
#' - DMR tool selection (combine multiple tools)
#' - Dynamic plots updating based on selections
#' - Quality control visualizations
#' - Downloadable results tables
#'
#' The report is rendered with Shiny runtime, allowing real-time updates
#' without recompilation.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Generate interactive report
#' methylkey_report(
#'   se = my_summarized_experiment,
#'   mrs = my_methylresultset,
#'   mvals_raw = raw_mvalues,
#'   mvals_sva = sva_corrected_mvalues,
#'   genome = "hg38",
#'   dmrtools = c("dmrcate", "dmrff"),
#'   output_file = "my_analysis.qmd",
#'   title = "Methylation Analysis Report",
#'   author = "Your Name"
#' )
#'
#' # Then render the report with:
#' # quarto::quarto_render("my_analysis.qmd")
#' }
#'
methylkey_report <- function(
  se,
  mrs,
  mvals_raw,
  mvals_sva = NULL,
  genome = "hg38",
  dmrtools = c("dmrcate", "dmrff"),
  output_file = "methylkey_report.qmd",
  title = "Methylation Analysis Report",
  author = "methylkey"
) {

  assertthat::assert_that(
    methods::is(mrs, "MethylResultSet"),
    msg = "mrs must be a MethylResultSet object"
  )
  assertthat::assert_that(
    is.matrix(mvals_raw),
    msg = "mvals_raw must be a matrix"
  )

  # Get metadata from MethylResultSet
  model <- mrs@metadata$model
  method <- mrs@metadata$method
  intercept <- mrs@metadata$intercept
  dmrtools_str <- paste(sprintf("'%s'", dmrtools), collapse = ", ")

  # Data filenames for the report to load
  mrs_file <- stringr::str_replace(output_file, ".qmd$", ".mrs.rds")
  mvals_raw_file <- stringr::str_replace(output_file, ".qmd$", ".mvals_raw.rds")
  mvals_sva_file <- stringr::str_replace(output_file, ".qmd$", ".mvals_sva.rds")
  se_file <- stringr::str_replace(output_file, ".qmd$", ".se.rds")

  # Build report content as character vector
  report_lines <- c(
    "---",
    paste0("title: \"", title, "\""),
    paste0("author: \"", author, "\""),
    "date: \"`r Sys.Date()`\"",
    "format:",
    "  html:",
    "    code-fold: true",
    "    code-summary: \"Show the code\"",
    "    code-tools: true",
    "    toc: true",
    "    toc-location: left",
    "    page-layout: full",
    "    lightbox: true",
    "    self-contained: true",
    "runtime: shiny",
    "server-type: default",
    "execute:",
    "  warning: false",
    "  error: false",
    "---",
    "",
    "```{r setup}",
    "#| include: false",
    "library(methylkey)",
    "library(tidyverse)",
    "library(DT)",
    "library(ggrepel)",
    "library(shiny)",
    "",
    "# Load pre-computed data",
    paste0("mrs <- readRDS('", mrs_file, "')"),
    paste0("mvals_raw <- readRDS('", mvals_raw_file, "')"),
    paste0("mvals_sva <- readRDS('", mvals_sva_file, "')"),
    paste0("se <- readRDS('", se_file, "')"),
    "betas <- methylkey::getBetas(se)",
    "pdata <- colData(se) %>% as.data.frame()",
    "",
    "params <- list(",
    paste0("  genome = '", genome, "',"),
    paste0("  dmrtools = c(", dmrtools_str, "),"),
    paste0("  model = '", model, "',"),
    paste0("  method = '", method, "',"),
    paste0("  intercept = '", intercept, "'"),
    ")",
    "```",
    "",
    "## Interactive Controls",
    "",
    "::: {layout-ncol=4}",
    "",
    "```{r model_selector}",
    "selectInput(",
    "  'selected_model',",
    "  'Select Model Contrast:',",
    "  choices = names(mrs@dmps),",
    "  selected = names(mrs@dmps)[1]",
    ")",
    "```",
    "",
    "```{r method_selector}",
    "selectInput(",
    "  'selected_method',",
    "  'Statistical Method:',",
    "  choices = c('ls', 'robust'),",
    paste0("  selected = '", method, "'"),
    ")",
    "```",
    "",
    "```{r data_selector}",
    "selectInput(",
    "  'data_source',",
    "  'Data Source:',",
    "  choices = c('Raw', 'SVA-corrected'),",
    "  selected = 'Raw'",
    ")",
    "```",
    "",
    "```{r tools_selector}",
    "checkboxGroupInput(",
    "  'selected_tools',",
    "  'DMR Tools:',",
    "  choices = params$dmrtools,",
    "  selected = params$dmrtools",
    ")",
    "```",
    "",
    ":::",
    "",
    "---",
    "",
    "## Quality Control",
    "",
    "### Sample Overview",
    "",
    "```{r sample_table}",
    "pdata %>%",
    "  DT::datatable(",
    "    extensions = 'Buttons',",
    "    options = list(",
    "      dom = 'Blfrtip',",
    "      buttons = c('copy', 'csv', 'excel')",
    "    )",
    "  )",
    "```",
    "",
    "### Distribution of Beta Values",
    "",
    "```{r beta_distribution}",
    "renderPlot({",
    "  betas_melted <- reshape2::melt(betas)",
    "  ggplot(betas_melted, aes(x = value, fill = Var2)) +",
    "    geom_density(alpha = 0.5) +",
    "    theme_minimal() +",
    "    labs(",
    "      x = 'Beta Value',",
    "      y = 'Density',",
    "      title = 'Distribution of Beta Values Across Samples'",
    "    ) +",
    "    theme(legend.position = 'none')",
    "})",
    "```",
    "",
    "---",
    "",
    "## Differential Methylation Analysis",
    "",
    "```{r dmp_analysis}",
    "renderUI({",
    "  index <- input$selected_model",
    "",
    "  # Get DMPs for selected model",
    "  dmps <- getDMPs(mrs, index)",
    "  sig_dmps <- dmps %>% filter(adj.P.Val < 0.05)",
    "",
    "  # Determine data source",
    "  if (input$data_source == 'SVA-corrected' && !is.null(mvals_sva)) {",
    "    mvals_use <- mvals_sva",
    "  } else {",
    "    mvals_use <- mvals_raw",
    "  }",
    "",
    "  tagList(",
    "    h3(paste('DMPs for:', index)),",
    "    p(paste('Number of significant DMPs (FDR < 0.05):', nrow(sig_dmps))),",
    "",
    "    tabsetPanel(",
    "      tabPanel(",
    "        'Table',",
    "        DT::datatable(",
    "          sig_dmps %>%",
    "            arrange(adj.P.Val) %>%",
    "            head(100),",
    "          extensions = 'Buttons',",
    "          options = list(",
    "            dom = 'Blfrtip',",
    "            buttons = c('copy', 'csv', 'excel')",
    "          )",
    "        )",
    "      ),",
    "      tabPanel(",
    "        'QQ Plot',",
    "        renderPlot({",
    "          qqman::qq(dmps$P.Value, main = paste('QQ plot:', index))",
    "        })",
    "      ),",
    "      tabPanel(",
    "        'Volcano Plot',",
    "        renderPlot({",
    "          my_volcano(dmps)",
    "        })",
    "      )",
    "    )",
    "  )",
    "})",
    "```",
    "",
    "### DMR Analysis",
    "",
    "```{r dmr_analysis}",
    "renderUI({",
    "  index <- input$selected_model",
    "  tools_use <- input$selected_tools",
    "",
    "  if (length(tools_use) == 0) {",
    "    return(p('Select at least one DMR tool', class = 'alert alert-warning'))",
    "  }",
    "",
    "  # Get results for selected tools and model",
    "  result <- getResults(",
    "    mrs,",
    "    index,",
    "    tools = tools_use,",
    "    genome = params$genome,",
    "    mvals = if (input$data_source == 'SVA-corrected') mvals_sva else mvals_raw",
    "  )",
    "",
    "  if (is.null(result) || nrow(result) == 0) {",
    "    return(p('No DMRs found with selected criteria'))",
    "  }",
    "",
    "  sig_dmrs <- result %>%",
    "    filter(fdr < 0.05, no.cpgs > 1) %>%",
    "    group_by(dmrtool, ID, fdr, genesUniq) %>%",
    "    summarise(n = n(), .groups = 'drop')",
    "",
    "  tagList(",
    "    h3(paste('DMRs for:', index)),",
    "    p(paste('Number of significant DMRs:', nrow(sig_dmrs))),",
    "",
    "    tabsetPanel(",
    "      tabPanel(",
    "        'Summary',",
    "        renderTable({",
    "          sig_dmrs %>%",
    "            group_by(dmrtool) %>%",
    "            summarise(",
    "              n_regions = n(),",
    "              .groups = 'drop'",
    "            )",
    "        })",
    "      ),",
    "      tabPanel(",
    "        'Table',",
    "        DT::datatable(",
    "          sig_dmrs %>%",
    "            arrange(fdr) %>%",
    "            head(100),",
    "          extensions = 'Buttons',",
    "          options = list(",
    "            dom = 'Blfrtip',",
    "            buttons = c('copy', 'csv', 'excel')",
    "          )",
    "        )",
    "      )",
    "    )",
    "  )",
    "})",
    "```"
  )

  # Write .qmd file
  writeLines(report_lines, output_file)

  # Save R data objects as RDS for the report to load
  saveRDS(mrs, mrs_file)
  saveRDS(mvals_raw, mvals_raw_file)
  if (!is.null(mvals_sva)) {
    saveRDS(mvals_sva, mvals_sva_file)
  }
  saveRDS(se, se_file)

  message(paste(
    "Interactive report created:",
    output_file,
    "\nRender with: quarto::quarto_render('", output_file, "')",
    sep = ""
  ))

  invisible(output_file)
}
