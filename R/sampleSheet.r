#' Format and rename columns in a sample sheet data frame
#'
#' This function formats and renames columns in a sample sheet
#'   data frame to ensure consistency and compatibility with
#'   downstream analysis. It performs the following operations:
#' - Converts column names to lowercase.
#' - Renames the "samples" or "sample_id" column to "samples".
#' - Renames the "barcode" or "basename" column to "barcode".
#' - Ensures that the "barcode" and "samples" columns exist in the sample sheet.
#' - Converts column names to valid R variable names.
#' - Converts numeric columns with discrete values to factors
#'   (integers with unique values < sqrt(n)).
#' - Converts character columns to factors.
#' - Separates the "barcode" column into "sentrix_id" and
#'   "sentrix_position" columns.
#' - Arranges the rows by the "barcode" column.
#' - Validates that "samples" and "barcode" columns contain unique values.
#'
#' @param sampleSheet A data frame representing the sample sheet.
#'
#' @return A processed and formatted sample sheet data frame.
#'
#' @importFrom assertthat assert_that
#' @importFrom dplyr rename rename_with everything mutate where across arrange
#' @importFrom tidyr separate
#'
#' @export
format_sample_sheet <- function(sample_sheet) {

  assertthat::assert_that(
    is.data.frame(sample_sheet),
    msg = "sampleSheet must be a data frame"
  )
  assertthat::assert_that(
    nrow(sample_sheet) > 0,
    msg = "sampleSheet cannot be empty"
  )

  # Convert column names to lowercase
  sample_sheet <- sample_sheet |>
    dplyr::rename_with(tolower)

  # Rename basename to barcode if present
  if ("basename" %in% colnames(sample_sheet)) {
    sample_sheet <- sample_sheet |>
      dplyr::rename(barcode = basename)
  }

  # Check that barcode column exists
  assertthat::assert_that(
    "barcode" %in% colnames(sample_sheet),
    msg = "barcode column is required (or basename as alternative)"
  )

  # Convert column names to valid R names and replace dots with underscores
  sample_sheet <- sample_sheet |>
    dplyr::rename_with(make.names) |>
    dplyr::rename_with(~gsub("\\.", "_", .), dplyr::everything())

  # Ensure first column is samples (rename if needed)
  first_col_name <- colnames(sample_sheet)[1]
  if (first_col_name != "samples") {
    sample_sheet <- sample_sheet |>
      dplyr::rename(samples = dplyr::any_of(first_col_name))
  }

  # Separate barcode into sentrix_id and sentrix_position
  if (!("sentrix_id" %in% colnames(sample_sheet)) &&
      !("sentrix_position" %in% colnames(sample_sheet))
  ) {
    sample_sheet <- sample_sheet |>
      tidyr::separate(
        .data$barcode,
        into = c("sentrix_id", "sentrix_position"),
        remove = FALSE,
        sep = "_",
        extra = "merge"
      )
  }

  # Convert numeric variables to factors if they are discrete
  # (all integers and number of unique values < sqrt(n))
  sample_sheet <- sample_sheet |>
    dplyr::mutate(dplyr::across(
      dplyr::where(is.numeric),
      ~if (all(. == as.integer(.), na.rm = TRUE) &&
          length(unique(na.omit(.))) < sqrt(length(.))
      ) {
        as.factor(.)
      } else {
        .
      }
    ))

  # Convert character columns to factors
  sample_sheet <- sample_sheet |>
    dplyr::mutate(dplyr::across(dplyr::where(is.character), as.factor))

  # Sort by barcode
  sample_sheet <- sample_sheet |>
    dplyr::arrange(.data$barcode)

  # Check that samples and barcode columns contain unique values
  assertthat::assert_that(
    length(unique(sample_sheet$samples)) == nrow(sample_sheet),
    msg = "samples column must contain unique values (no duplicates)"
  )
  assertthat::assert_that(
    length(unique(sample_sheet$barcode)) == nrow(sample_sheet),
    msg = "barcode column must contain unique values (no duplicates)"
  )

  sample_sheet
}