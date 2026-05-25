# Extracted from test_newBetas.R:16

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "methylkey", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
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
