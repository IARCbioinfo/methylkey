test_that("format_sample_sheet handles basic input correctly", {
  # Create a simple sample sheet
  ss <- data.frame(
    SampleID = c("Sample1", "Sample2", "Sample3"),
    Barcode = c("203021070069_R03C01", "203021070069_R04C01", "203021070069_R05C01"),
    Group = c("Control", "Case", "Control"),
    stringsAsFactors = FALSE
  )
  
  result <- format_sample_sheet(ss)
  
  # Check that result is a data frame
  expect_true(is.data.frame(result))
  
  # Check that required columns exist
  expect_true("barcode" %in% colnames(result))
  expect_true("samples" %in% colnames(result))
  
  # Check that sentrix columns are created
  expect_true("sentrix_id" %in% colnames(result))
  expect_true("sentrix_position" %in% colnames(result))
  
  # Check that samples are sorted by barcode
  expect_equal(result$barcode, sort(result$barcode))
})

test_that("format_sample_sheet renames basename to barcode", {
  ss <- data.frame(
    SampleID = c("Sample1", "Sample2"),
    Basename = c("203021070069_R03C01", "203021070069_R04C01"),
    stringsAsFactors = FALSE
  )
  
  result <- format_sample_sheet(ss)
  expect_true("barcode" %in% colnames(result))
  expect_false("basename" %in% colnames(result))
})

test_that("format_sample_sheet rejects empty data frame", {
  ss <- data.frame()
  expect_error(format_sample_sheet(ss), "cannot be empty")
})

test_that("format_sample_sheet requires barcode column", {
  ss <- data.frame(
    SampleID = c("Sample1", "Sample2"),
    stringsAsFactors = FALSE
  )
  expect_error(format_sample_sheet(ss), "barcode column is required")
})
