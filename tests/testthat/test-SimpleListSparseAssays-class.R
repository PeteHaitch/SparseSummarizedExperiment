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

test_that("SimpleListSparseAssays class definition is correct", {
  expect_true(extends("SimpleListSparseAssays", "SparseAssays"))
  expect_true(extends("SimpleListSparseAssays", "SimpleList"))
})

test_that(".valid.SimpleListSparseAssays and validObject work", {
  expect_null(.valid.SimpleListSparseAssays(SimpleList()))

  msg <- paste0("All sample-level data within each sparse assay of a ",
                "'SimpleListSparseAssays' object must have one element ",
                "named 'key', one element named 'value', and nothing ",
                "else.")
  bad_sa <- sa
  bad_sa@listData$sa1$s1 <- SimpleList(dog = 1:10, cat = matrix(1:20, ncol = 2))
  expect_equal(.valid.SimpleListSparseAssays(bad_sa), msg)
  expect_error(validObject(bad_sa), msg)

  msg <- paste0("All 'key' elements of a 'SimpleListSparseAssays' object must ",
                "be integer vectors.")
  bad_sa <- sa
  bad_sa@listData$sa1$s1$key <- as.numeric(1:10)
  expect_identical(.valid.SimpleListSparseAssays(bad_sa), msg)
  expect_error(validObject(bad_sa), msg)

  msg <- paste0("All 'value' elements of a 'SimpleListSparseAssays' object ",
                "must be numeric matrix objects.")
  bad_sa <- sa
  bad_sa@listData$sa1$s1$value <- as.data.frame(bad_sa@listData$sa1$s1$value)
  expect_identical(.valid.SimpleListSparseAssays(bad_sa), msg)
  expect_error(validObject(bad_sa), msg)
  bad_sa <- sa
  bad_sa@listData$sa1$s1$value <- matrix("a", nrow = nrow(bad_sa), ncol = 2)
  expect_identical(.valid.SimpleListSparseAssays(bad_sa), msg)
  expect_error(validObject(bad_sa), msg)

  msg <- paste0("All sparse assays of a 'SimpleListSparseAssays' object must ",
                "have an identical number of samples.")
  bad_sa <- c(sa, sa[, 2])
  expect_identical(.valid.SimpleListSparseAssays(bad_sa), msg)
  expect_error(validObject(bad_sa), msg)

  msg <- paste0("All sparse assays of a 'SimpleListSparseAssays' object must ",
                "have identical sample names.")
  bad_sa <- c(sa, sa)
  expect_error(names(bad_sa[[1]]) <- paste0("S", 1:6),
               msg)
  names(bad_sa@listData[[1]]) <- paste0("S", 1:6)
  expect_identical(.valid.SimpleListSparseAssays(bad_sa), msg)
  expect_error(validObject(bad_sa), msg)

  msg <- paste0("All 'key' elements of a 'SimpleListSparseAssays' object must ",
                "have identical length.")
  bad_sa <- sa
  bad_sa@listData$sa1$s1$key <- 1:9
  expect_identical(.valid.SimpleListSparseAssays(bad_sa), msg)
  expect_error(validObject(bad_sa), msg)
  bad_sa <- c(sa, sa[1:9, ])
  expect_identical(.valid.SimpleListSparseAssays(bad_sa), msg)
  expect_error(validObject(bad_sa), msg)

  msg <- paste0("All 'data' elements within each sparse assay of a ",
                "'SimpleListSparseAssays' object must have identical ncol.")
  bad_sa <- sa
  bad_sa@listData$sa1$s1$value <- matrix(1:5, ncol = 1)
  expect_identical(.valid.SimpleListSparseAssays(bad_sa), msg)
  expect_error(validObject(bad_sa), msg)

  msg <- paste0("Maximum value in each 'key' element must be less than or ",
         "equal to the number of rows in each corresponding 'value' ",
         "element of a 'SimpleListSparseAssays' object.")
  bad_sa <- sa
  bad_sa@listData$sa1$s1$key <- rep(9999999L, length(bad_sa@listData$sa1$s1$key))
  expect_identical(.valid.SimpleListSparseAssays(bad_sa), msg)
  expect_error(validObject(bad_sa), msg)

  expect_null(.valid.SimpleListSparseAssays(sa))
  expect_true(validObject(sa))
})

test_that("dim,SimpleListSparseAssays-method works", {
  expect_identical(dim(SparseAssays()), c(0L, 0L))
  expect_identical(dim(sa), c(10000L, 6L))
  expect_identical(dim(sa[1:10, 3]), c(10L, 1L))
})

test_that("[,SimplelistSparseAssays-method works", {

  x1 <- SimpleList(
    s1 = SimpleList(key = as.integer(c(NA, 1, NA, NA, 2, NA, 3, NA, 4, 5)),
                    value = matrix(1:10, ncol = 2)),
    s2 = SimpleList(key = as.integer(c(NA, NA, 1, 2, NA, NA, 3, 4, NA, NA)),
                    value = matrix(8:1, ncol = 2)))
  y1 <- SimpleList(
    s1 = SimpleList(key = as.integer(c(NA, 1, NA, NA, 2, NA, 3)),
                    value = matrix(c(1:3, 6:8), ncol = 2)),
    s2 = SimpleList(key = as.integer(c(NA, NA, 1, 2, NA, NA, 3)),
                    value = matrix(c(8:6, 4:2), ncol = 2)))
  z1 <- SimpleList(
    s2 = SimpleList(key = as.integer(c(NA, NA, 1, 2, NA, NA, 3, 4, NA, NA)),
                    value = matrix(8:1, ncol = 2)))
  w1 <- SimpleList(
    s2 = SimpleList(key = as.integer(c(NA, NA, 1, 2, NA, NA, 3)),
                    value = matrix(c(8:6, 4:2), ncol = 2)))
  x2 <- SimpleList(
    s1 = SimpleList(key = as.integer(c(NA, 1, NA, 2, 2, NA, 1, NA, NA, 1)),
                    value = matrix(1:2, ncol = 1)),
    s2 = SimpleList(key = as.integer(c(1, 1, 1, 2, NA, NA, NA, NA, NA, NA)),
                    value = matrix(4:3, ncol = 1)))
  y2 <- SimpleList(
    s1 = SimpleList(key = as.integer(c(NA, 1, NA, 2, 2, NA, 1)),
                    value = matrix(1:2, ncol = 1)),
    s2 = SimpleList(key = as.integer(c(1, 1, 1, 2, NA, NA, NA)),
                    value = matrix(4:3, ncol = 1)))
  z2 <- SimpleList(
    s2 = SimpleList(key = as.integer(c(1, 1, 1, 2, NA, NA, NA, NA, NA, NA)),
                    value = matrix(4:3, ncol = 1)))
  w2 <- SimpleList(
    s2 = SimpleList(key = as.integer(c(1, 1, 1, 2, NA, NA, NA)),
                    value = matrix(4:3, ncol = 1)))
  x <- SparseAssays(SimpleList(sa1 = x1, sa2 = x2))
  y <- SparseAssays(SimpleList(sa1 = y1, sa2 = y2))
  z <- SparseAssays(SimpleList(sa1 = z1, sa2 = z2))
  w <- SparseAssays(SimpleList(sa1 = w1, sa2 = w2))

  # NOTE: Whenever [,SimpleListSparseAssays method is called with 'i', the
  #       value slot will be sorted using data.table::setkey(), whose
  #       sort-order is equivalent to base::sort(x, na.last = FALSE). There are
  #       multiple, equivalent definitions of SimpleListSparseAssays data.
  #       These are equivalent up to rows permutations.
  expect_that(x, not(is_identical_to(x[seq_len(nrow(x)), ])))
  expect_identical(densify(x, seq_along(x), seq_len(ncol(x))),
                   densify(x[seq_len(nrow(x)), ],
                           seq_along(x[seq_len(nrow(x)), ]),
                           seq_len(ncol(x[seq_len(nrow(x)), ]))))
  expect_that(x[1:7, ], not(is_identical_to(y)))
  expect_identical(densify(x[1:7, ],
                           seq_along(x[1:7, ]),
                           seq_len(ncol(x[1:7, ]))),
                   densify(y,
                           seq_along(y[seq_len(nrow(y)), ]),
                           seq_len(ncol(y[seq_len(nrow(y)), ]))))
  # NOTE: Should be identical when only j is specified
  expect_that(x[, 2], is_identical_to(z))
  expect_identical(densify(x[, 2],
                           seq_along(x[, 2]),
                           seq_len(ncol(x[, 2]))),
                   densify(z,
                           seq_along(z[seq_len(nrow(z)), ]),
                           seq_len(ncol(z[seq_len(nrow(z)), ]))))
  expect_that(x[1:7, 2], not(is_identical_to(w)))
  expect_identical(densify(x[1:7, 2],
                           seq_along(x[1:7, 2]),
                           seq_len(ncol(x[1:7, 2]))),
                   densify(w,
                           seq_along(w[seq_len(nrow(w)), ]),
                           seq_len(ncol(w[seq_len(nrow(w)), ]))))
  # Testing non-contiguous i and j.
  x <- matrix(1:10, ncol = 2, dimnames = list(letters[1:5], LETTERS[1:2]))
  xx <- SparseAssays(x)
  names(xx[[1L]]) <- "X"
  y <- matrix(101:110, ncol = 2, dimnames = list(letters[1:5], LETTERS[1:2]))
  yy <- SparseAssays(y)
  names(yy[[1L]]) <- "Y"
  z <- matrix(1001:1010, ncol = 2, dimnames = list(letters[1:5], LETTERS[1:2]))
  zz <- SparseAssays(z)
  names(zz[[1L]]) <- "Z"
  w <- cbind(xx, yy, zz)
  expect_identical(densify(xx[c(1, 4, 2), ], 1, 1)[[1]][[1]],
                   matrix(as.integer(c(1, 4, 2, 6, 9, 7)), ncol = 2,
                          dimnames = list(c("a", "d", "b"), c("A", "B"))))
  expect_identical(w[, c(3, 1)], cbind(zz, xx))
  expect_identical(densify(w[c(1, 4, 2), c(3, 1)], 1, 1:2),
                   SimpleList(
                     list(list(Z = matrix(as.integer(c(1, 4, 2, 6, 9, 7)) + 1000L,
                                          ncol = 2,
                                          dimnames = list(c("a", "d", "b"),
                                                          c("A", "B"))),
                               X = matrix(as.integer(c(1, 4, 2, 6, 9, 7)),
                                          ncol = 2,
                                          dimnames = list(c("a", "d", "b"),
                                                          c("A", "B")))))))

  # Testing non-numeric i and j
  expect_identical(densify(xx[c("a", "d", "b"), ], 1, 1)[[1]][[1]],
                   matrix(as.integer(c(1, 4, 2, 6, 9, 7)), ncol = 2,
                          dimnames = list(c("a", "d", "b"), c("A", "B"))))
  expect_identical(w[, c(3, 1)], cbind(zz, xx))
  expect_identical(densify(w[c("a", "d", "b"), c("Z", "X")], 1, 1:2),
                   SimpleList(
                     list(list(Z = matrix(as.integer(c(1, 4, 2, 6, 9, 7)) + 1000L,
                                          ncol = 2,
                                          dimnames = list(c("a", "d", "b"),
                                                          c("A", "B"))),
                               X = matrix(as.integer(c(1, 4, 2, 6, 9, 7)),
                                          ncol = 2,
                                          dimnames = list(c("a", "d", "b"),
                                                          c("A", "B")))))))
  # Test out-of-bounds indices
  msg <- "subscript contains NAs or out-of-bounds indices"
  expect_error(xx[nrow(xx) + 1, ], msg)
  expect_error(xx["bad_idx", ], msg)
  expect_error(xx[nrow(x) + 1, 1], msg)
  expect_error(xx["bad_idx", 1], msg)
  expect_error(xx[, ncol(x) + 1], msg)
  expect_error(xx[, "bad_idx"], "subscript contains invalid names")
  expect_error(xx[nrow(x) + 1, ncol(x) + 1], msg)
  expect_error(xx["bad_idx", "another_bad_idx"], msg)

  # Trigger warning re drop argument
  expect_warning(xx[1, , drop = TRUE],
                 "'drop' ignored '\\[,SimpleListSparseAssays,ANY,ANY-method'")
})

# Test all code in R/SimpleListSparseAssays-class.R
test_that("UP TO HERE: Test [<-, then rest of R/SimpleListSparseAssays-class.R", {
  expect_true(TRUE)
})
