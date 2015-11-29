context("Functionality that probably belongs in the GenomicRanges package")

test_that("combine,GRanges,GRanges-method is compatible with generic", {
  generic <- getGeneric("combine")
  method <- getMethod("combine", c("GRanges", "GRanges"))
  expect_identical(generic@signature, c("x", "y"))
  expect_identical(formals(generic@.Data), formals(method@.Data))
})

test_that("combine,GRanges,GRanges-method works", {
  set.seed(666)
  x <- simGR(10)
  expect_identical(combine(x, GRanges()), x)
  expect_identical(combine(GRanges(), x), x)
  expect_identical(combine(x, x, x), x)
  expect_identical(combine(x[1:7], x[3:9]), x[1:9])
  expect_identical(combine(x[10:5], x[4:1]), x[10:1])
})

test_that("combine,GRangesList,GRangesList-method is compatible with generic", {
  generic <- getGeneric("combine")
  method <- getMethod("combine", c("GRangesList", "GRangesList"))
  expect_identical(generic@signature, c("x", "y"))
  expect_identical(formals(generic@.Data), formals(method@.Data))
})


test_that("combine,GRangesList,GRangesList-method works", {
  set.seed(666)
  gr <- simGR(10)
  x <- GenomicRanges::GRangesList(gr[1:4], gr[5:10])
  expect_identical(combine(x, GRangesList()), x)
  expect_identical(combine(GRangesList(), x), x)
  x1 <- x[1]
  x2 <- x[2]
  expect_error(combine(x1, x2),
               paste0("'names' of 'x' and 'y' must be non-NULL when combining ",
                      "'GRangesList' object"))
  xn <- x
  names(xn) <- c("A", "B")
  x1 <- xn[1]
  x2 <- xn[2]
  expect_identical(combine(x1, x2), xn)
  expect_identical(combine(x1, x1), x1)
  x1a <- GenomicRanges::GRangesList(A = gr[1:3])
  x1b <- GenomicRanges::GRangesList(A = gr[2:4])
  expect_identical(combine(x1a, x1b), x1)
})
