test_that("CpGNAexcl identifies probes with many NAs", {
  # Create a test matrix with some NAs
  betas <- matrix(runif(100), nrow = 20, ncol = 5)
  rownames(betas) <- paste0("cg", sprintf("%07d", 1:20))
  
  # Add NAs to first row (30% NAs)
  betas[1, 1:2] <- NA
  
  # Add NAs to second row (40% NAs)
  betas[2, 1:2] <- NA
  
  # Test with nalimit = 0.2 (20%)
  excluded <- CpGNAexcl(betas, nalimit = 0.2)
  
  # First two rows have >= 20% NAs, should be excluded
  expect_true(rownames(betas)[1] %in% excluded)
  expect_true(rownames(betas)[2] %in% excluded)
  
  # Other rows should not be excluded
  expect_false(rownames(betas)[3] %in% excluded)
})

test_that("CpGNAexcl handles matrix with no NAs", {
  betas <- matrix(runif(100), nrow = 20, ncol = 5)
  rownames(betas) <- paste0("cg", sprintf("%07d", 1:20))
  
  excluded <- CpGNAexcl(betas, nalimit = 0.2)
  
  # No probes should be excluded
  expect_equal(length(excluded), 0)
})

test_that("CpGNAexcl respects nalimit parameter", {
  betas <- matrix(runif(100), nrow = 10, ncol = 5)
  rownames(betas) <- paste0("cg", sprintf("%06d", 1:10))
  
  # Add 40% NAs to first row
  betas[1, 1:2] <- NA
  
  # With nalimit = 0.5, first row should NOT be excluded
  excluded_50 <- CpGNAexcl(betas, nalimit = 0.5)
  expect_false(rownames(betas)[1] %in% excluded_50)
  
  # With nalimit = 0.2, first row should be excluded
  excluded_20 <- CpGNAexcl(betas, nalimit = 0.2)
  expect_true(rownames(betas)[1] %in% excluded_20)
})
