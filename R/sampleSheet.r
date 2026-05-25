#' Format and rename columns in a sample sheet data frame
#'
#' This function formats and renames columns in a sample sheet data frame to ensure consistency and compatibility with downstream analysis. It performs the following operations:
#' - Converts column names to lowercase.
#' - Renames the "samples" or "sample_id" column to "samples".
#' - Renames the "barcode" or "basename" column to "barcode".
#' - Ensures that the "barcode" and "samples" columns exist in the sample sheet.
#' - Converts column names to valid R variable names.
#' - Converts numeric columns with discrete values to factors (integers with unique values < sqrt(n)).
#' - Converts character columns to factors.
#' - Separates the "barcode" column into "sentrix_id" and "sentrix_position" columns.
#' - Arranges the rows by the "barcode" column.
#' - Validates that "samples" and "barcode" columns contain unique values.
#'
#' @param sampleSheet A data frame representing the sample sheet.
#'
#' @return A processed and formatted sample sheet data frame.
#'
#'
#' @export
formatSampleSheet <- function(sampleSheet) {

  assertthat::assert_that(
    is.data.frame(sampleSheet),
    msg = "sampleSheet must be a data frame"
  )
  assertthat::assert_that(
    nrow(sampleSheet) > 0,
    msg = "sampleSheet cannot be empty"
  )

  # Convert column names to lowercase
  sampleSheet <- sampleSheet |>
    dplyr::rename_with(tolower)

  # Rename basename to barcode if present
  if ("basename" %in% colnames(sampleSheet)) {
    sampleSheet <- sampleSheet |>
      dplyr::rename(barcode = basename)
  }

  # Check that barcode column exists
  assertthat::assert_that(
    "barcode" %in% colnames(sampleSheet),
    msg = "barcode column is required (or basename as alternative)"
  )

  # Convert column names to valid R names and replace dots with underscores
  sampleSheet <- sampleSheet |>
    dplyr::rename_with(make.names) |>
    dplyr::rename_with(~gsub("\\.", "_", .), dplyr::everything())

  # Ensure first column is samples (rename if needed)
  first_col_name <- colnames(sampleSheet)[1]
  if (first_col_name != "samples") {
    sampleSheet <- sampleSheet |>
      dplyr::rename(samples = dplyr::all_of(first_col_name))
  }

  # Separate barcode into sentrix_id and sentrix_position if not already separated
  if (!("sentrix_id" %in% colnames(sampleSheet)) &&
      !("sentrix_position" %in% colnames(sampleSheet))) {
    sampleSheet <- sampleSheet |>
      tidyr::separate(
        barcode,
        into = c("sentrix_id", "sentrix_position"),
        remove = FALSE,
        sep = "_",
        extra = "merge"
      )
  }

  # Convert numeric variables to factors if they are discrete
  # (all integers and number of unique values < sqrt(n))
  sampleSheet <- sampleSheet |>
    dplyr::mutate(dplyr::across(
      dplyr::where(is.numeric),
      ~if (all(. == as.integer(.), na.rm = TRUE) &&
            length(unique(na.omit(.))) < sqrt(length(.))) {
        as.factor(.)
      } else {
        .
      }
    ))

  # Convert character columns to factors
  sampleSheet <- sampleSheet |>
    dplyr::mutate(dplyr::across(dplyr::where(is.character), as.factor))

  # Sort by barcode
  sampleSheet <- sampleSheet |>
    dplyr::arrange(barcode)

  # Check that samples and barcode columns contain unique values
  assertthat::assert_that(
    length(unique(sampleSheet$samples)) == nrow(sampleSheet),
    msg = "samples column must contain unique values (no duplicates)"
  )
  assertthat::assert_that(
    length(unique(sampleSheet$barcode)) == nrow(sampleSheet),
    msg = "barcode column must contain unique values (no duplicates)"
  )

  return(sampleSheet)
}