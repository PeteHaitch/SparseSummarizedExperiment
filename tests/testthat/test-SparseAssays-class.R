context("SparseAssays class")

test_that("NROW,SimpleListSparseAssays-method is compatible with generic", {
  generic <- getGeneric("NROW")
  method <- getMethod("NROW", "SparseAssays")
  expect_identical(generic@signature, "x")
  expect_identical(formals(generic@.Data), formals(method@.Data))
})
