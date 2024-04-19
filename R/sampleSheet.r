#' Format and rename columns in a sample sheet data frame
#'
#' This function formats and renames columns in a sample sheet data frame to ensure consistency and compatibility with downstream analysis. It performs the following operations:
#' - Converts column names to lowercase.
#' - Renames the "samples" or "sample_id" column to "samples".
#' - Renames the "barcode" or "basename" column to "barcode".
#' - Ensures that the "barcode" and "samples" columns exist in the sample sheet.
#' - Converts column names to valid R variable names.
#' - Converts character columns to factors.
#' - Separates the "barcode" column into "sentrix_id" and "sentrix_position" columns.
#' - Arranges the rows by the "barcode" column.
#'
#' @param sampleSheet A data frame representing the sample sheet.
#'
#' @return A processed and formatted sample sheet data frame.
#'
#'
#' @export
formatSampleSheet<-function(sampleSheet){
  
  sampleSheet<-sampleSheet %>% rename_with(tolower) %>% 
  #  dplyr::rename(any_of(c(samples="sample",samples="sample_id"))) %>%
    dplyr::rename(any_of(c(barcode="basename")))
  assertthat::assert_that( "barcode" %in% colnames(sampleSheet), msg="The barcode column is missing !" )
  #assertthat::assert_that( "samples" %in% colnames(sampleSheet), msg="The samples column is missing !" )
  
  sampleSheet <- sampleSheet  %>% 
    rename_with( make.names, unique = TRUE ) %>%
    setNames(gsub("\\.", "_", names(.))) %>%
    dplyr::rename(samples=1) %>% # first column must be samples names
    tidyr::separate(barcode, into = c("sentrix_id","sentrix_position"), remove=F) %>%
    dplyr::mutate_if( is.character, as.factor ) %>%
    dplyr::arrange(barcode)
}