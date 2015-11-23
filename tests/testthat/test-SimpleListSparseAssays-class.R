context("SimpleListSparseAssays class")

test_that(".valid.SparseAssays() works", {
  expect_null(.valid.SparseAssays(sa))
  # NOTE: Not much to test here, almost anything it seems can be coerced to
  #       a SimpleList. It means that this validity check is rather useless
  #       (but recall that subclasses will/should have their own additional
  #       validity checks).
  expect_null(.valid.SparseAssays(matrix()))
  expect_null(.valid.SparseAssays(SimpleList()))
})

test_that(".normarg.sparse_assays() works", {
  expect_is(.normarg.sparse_assays(S4Vectors::SimpleList()), "SimpleList")
  expect_is(.normarg.sparse_assays(list()), "SimpleList")
  expect_is(.normarg.sparse_assays(matrix()), "SimpleList")
  sa_matrix <- .normarg.sparse_assays(matrix(1:10, ncol = 2))
  expect_identical(names(sa_matrix[[1L]][[1L]]), c("key", "value"))
  expect_error(.normarg.sparse_assays(1:10),
               "'sparse_assays' must be a SimpleList, list, or matrix")
})

test_that("SparseAssays() works", {
  sa_matrix <- SparseAssays(matrix(1:10, ncol = 2))
  expect_is(sa_matrix, "SparseAssays")
  expect_is(sa_matrix, "SimpleListSparseAssays")
  expect_error(SparseAssays(matrix(1:10, ncol = 2), subclass = "FancyClass"),
               "'subclass' error: 'FancyClass' does not extend 'SparseAssays'")
  expect_error(SparseAssays(matrix(1:10, ncol = 2),
                            subclass = "ShallowSimpleListAssays"),
               paste0("'subclass' error: 'ShallowSimpleListAssays' does not ",
                      "extend 'SparseAssays'"))
})

test_that("length,SimpleListSparseAssays-method works", {
  expect_identical(length(sa), 1L)
  expect_identical(length(c(sa, sa, sa, sa, sa)), 5L)
})

test_that("NROW,SimpleListSparesAssays-method works", {
  expect_identical(NROW(sa), dim(sa)[1L])
})

test_that("names,SimpleListSparseAssays-method works", {
  expect_identical(names(sa), "sa1")
})

test_that("names<-,SimpleListSparseAssays-method works", {
  SA <- sa
  names(SA) <- "SA1"
  expect_identical(names(SA), "SA1")
})

test_that("[[,SimpleListSparseAssays-method works", {
  sa1 <- sa[[1L]]
  expect_is(sa1, "SimpleList")
  expect_identical(names(sa1), paste0("s", 1:6))
  x <- c(sa, sa)
  names(x) <- c("A", "sa1")
  expect_identical(x[[2L]], sa1)
})

test_that("[[<-,SimpleListSparseAssays-method works", {
  x <- sa[1:100, ]
  y <- c(x, x)
  expect_that(y[[2L]] <- sa[101:200, ][[1L]], not(throws_error()))
  expect_is(y, "SparseAssays")
  expect_error(y[[2L]] <- sa[101:300, ][[1L]],
               paste0("All 'key' elements of a 'SimpleListSparseAssays' ",
                      "object must have identical length."))
  expect_error(y[[2L]] <- sa[101:200, ],
               paste0("All sample-level data within each sparse assay of a ",
                      "'SimpleListSparseAssays' object must have one element ",
                      "named 'key', one element named 'value', and nothing ",
                      "else."))
})

test_that("UP TO HERE", {
  # Test all code in R/SimpleListSparseAssays-class.R
  expect_true(FALSE)
})
