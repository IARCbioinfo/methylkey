test_that("getPlateform detects EPIC correctly", {
  # Create a matrix with EPIC number of probes
  mat_epic <- matrix(0, nrow = 866553, ncol = 5)
  result <- getPlateform(mat_epic)
  expect_equal(result, "IlluminaHumanMethylationEPIC")
})

test_that("getPlateform detects 450k correctly", {
  # Create a matrix with 450k number of probes
  mat_450k <- matrix(0, nrow = 486427, ncol = 5)
  result <- getPlateform(mat_450k)
  expect_equal(result, "IlluminaHumanMethylation450k")
})

test_that("getPlateform detects Mouse 285k correctly", {
  # Create a matrix with mouse array number of probes
  mat_mouse <- matrix(0, nrow = 296070, ncol = 5)
  result <- getPlateform(mat_mouse)
  expect_equal(result, "IlluminaMouseMethylation285k")
})

test_that("getPlateform returns unknown for unmapped size", {
  mat_unknown <- matrix(0, nrow = 999999, ncol = 5)
  result <- getPlateform(mat_unknown)
  expect_equal(result, "unknown")
})

test_that("getPlateform handles NULL input", {
  result <- getPlateform(NULL)
  expect_equal(result, "unknown")
})
