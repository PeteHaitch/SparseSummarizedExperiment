context("Functionality that probably belongs in the GenomicRanges package")

test_that("combine,GenomicRanges,GenomicRanges-method is compatible with generic", {
  generic <- getGeneric("combine")
  method <- getMethod("combine", c("GenomicRanges", "GenomicRanges"))
  expect_identical(generic@signature, c("x", "y"))
  expect_identical(formals(generic@.Data), formals(method@.Data))
})

test_that("combine,GenomicRanges,GenomicRanges-method works on duplicate-free inputs", {
  set.seed(666)
  x <- simGR(10)
  expect_identical(combine(x, GRanges()), x)
  expect_identical(combine(GRanges(), x), x)
  expect_identical(combine(x, x, x), x)
  expect_identical(combine(x[1:7], x[3:9]), x[1:9])
  expect_identical(combine(x[10:5], x[4:1]), x[10:1])
})

test_that("combine,GenomicRanges,GenomicRanges-method handles duplicate elements", {
  set.seed(666)
  x <- simGR(10)
  expect_identical(combine(x[c(1, 1)], x[1]), x[1])
})

test_that("combine,GenomicRanges,GenomicRanges-method errors on 'non-identical' duplicate elements", {
  set.seed(666)
  x <- simGR(10)
  xx <- x[c(1, 1)]
  mcols(xx)[2, 1] <- mcols(x)[2, 1]
  expect_error(combine(xx, x),
               "x contains duplicate ranges whose 'mcols\\(\\)' differ")
  expect_error(combine(x, xx),
               "y contains duplicate ranges whose 'mcols\\(\\)' differ")
  expect_error(suppressWarnings(combine(xx[1], xx[2])),
               "'mcols\\(x\\)' and 'mcols\\(y\\)' are not compatible.")
  expect_warning(try(combine(xx[1], xx[2]), silent = TRUE),
                 "data frame column 'feature_id' shared rows not all equal")
})
