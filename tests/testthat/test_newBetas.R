test_that("newBetas creates valid Betas object", {
  # Create simple test data
  betas_matrix <- matrix(runif(100), nrow = 20, ncol = 5)
  colnames(betas_matrix) <- c("Sample1", "Sample2", "Sample3", "Sample4", "Sample5")
  rownames(betas_matrix) <- paste0("cg", sprintf("%07d", 1:20))
  
  ss <- data.frame(
    sample_name = c("Sample1", "Sample2", "Sample3", "Sample4", "Sample5"),
    barcode = c("203021070069_R03C01", "203021070069_R04C01", 
                "203021070069_R05C01", "203021070069_R06C01",
                "203021070069_R07C01"),
    group = c("A", "A", "B", "B", "B"),
    stringsAsFactors = FALSE
  )
  
  result <- newBetas(betas_matrix, ss, na = 0.2)
  
  # Check class
  expect_true(class(result) == "Betas")
  
  # Check dimensions
  expect_equal(nrow(result), 20)
  expect_equal(ncol(result), 5)
  
  # Check that metadata is present
  expect_true(!is.null(metadata(result)$plateform))
})

test_that("newBetas rejects non-matrix input", {
  betas_df <- data.frame(matrix(runif(100), nrow = 20, ncol = 5))
  ss <- data.frame(
    sample_name = c("S1", "S2", "S3", "S4", "S5"),
    barcode = c("203021070069_R03C01", "203021070069_R04C01", 
                "203021070069_R05C01", "203021070069_R06C01",
                "203021070069_R07C01")
  )
  
  expect_error(newBetas(betas_df, ss, na = 0.2), "must be a matrix")
})

test_that("newBetas requires barcode column", {
  betas_matrix <- matrix(runif(100), nrow = 20, ncol = 5)
  colnames(betas_matrix) <- c("S1", "S2", "S3", "S4", "S5")
  
  ss <- data.frame(
    sample_name = c("S1", "S2", "S3", "S4", "S5"),
    stringsAsFactors = FALSE
  )
  
  expect_error(newBetas(betas_matrix, ss, na = 0.2), "barcode column")
})

test_that("newBetas validates na parameter", {
  betas_matrix <- matrix(runif(100), nrow = 20, ncol = 5)
  colnames(betas_matrix) <- c("S1", "S2", "S3", "S4", "S5")
  
  ss <- data.frame(
    sample_name = c("S1", "S2", "S3", "S4", "S5"),
    barcode = c("203021070069_R03C01", "203021070069_R04C01", 
                "203021070069_R05C01", "203021070069_R06C01",
                "203021070069_R07C01")
  )
  
  expect_error(newBetas(betas_matrix, ss, na = 1.5), "must be between 0 and 1")
  expect_error(newBetas(betas_matrix, ss, na = -0.1), "must be between 0 and 1")
})
