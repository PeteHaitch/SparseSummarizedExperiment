context("SimpleListSparseAssays class")

# NOTE: x, y, z, and w are defined in helper-make-test-data.R, as is
#       identical_SparseAssays().

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

test_that("[,SimplelistSparseAssays-method works (no index)", {
  # NOTE: Whenever [,SimpleListSparseAssays method is called with 'i', the
  #       value slot will be sorted using data.table::setkey(), whose
  #       sort-order is equivalent to base::sort(x, na.last = FALSE). There are
  #       multiple, equivalent definitions of SimpleListSparseAssays data.
  #       These are equivalent up to rows permutations.
  # NOTE: While this might properly be considered a bug, subsetting
  #       SparseAssays should be consistent with subsetting Assays in that at
  #       least one of i or j must be specified.
  expect_error(x[], "object 'fun' not found")
  expect_error(as(x, "ShallowSimpleListAssays")[], "object 'fun' not found")
})

test_that("[,SimplelistSparseAssays-method works (i only)", {
  expect_that(x, not(is_identical_to(x[seq_len(nrow(x)), ])))
  expect_true(identical_SparseAssays(x, x[seq_len(nrow(x)), ]))
  expect_that(x[1:7, ], not(is_identical_to(y)))
  expect_true(identical_SparseAssays(x[1:7, ], y))
})

test_that("[,SimplelistSparseAssays-method works (j only)", {
  expect_that(x[, 2], is_identical_to(z))
  expect_true(identical_SparseAssays(x[, 2], z))
})

test_that("[,SimplelistSparseAssays-method works (i and j)", {
  expect_that(x[1:7, 2], not(is_identical_to(w)))
  expect_true(identical_SparseAssays(x[1:7, 2], w))
})

test_that("[,SimplelistSparseAssays-method works (non-contiguous i and j)", {
  expect_identical(densify(X[c(1, 4, 2), ], 1, 1)[[1]][[1]],
                   matrix(as.integer(c(1, 4, 2, 6, 9, 7)), ncol = 2,
                          dimnames = list(c("a", "d", "b"), c("A", "B"))))
  expect_identical(W[, c(3, 1)], cbind(Z, X))
  expect_identical(densify(W[c(1, 4, 2), c(3, 1)], 1, 1:2),
                   SimpleList(
                     list(list(Z = matrix(as.integer(c(1, 4, 2, 6, 9, 7)) + 1000L,
                                          ncol = 2,
                                          dimnames = list(c("a", "d", "b"),
                                                          c("A", "B"))),
                               X = matrix(as.integer(c(1, 4, 2, 6, 9, 7)),
                                          ncol = 2,
                                          dimnames = list(c("a", "d", "b"),
                                                          c("A", "B")))))))
})

test_that("[,SimplelistSparseAssays-method works (non-contiguous i and j)", {
  expect_identical(densify(X[c("a", "d", "b"), ], 1, 1)[[1]][[1]],
                   matrix(as.integer(c(1, 4, 2, 6, 9, 7)), ncol = 2,
                          dimnames = list(c("a", "d", "b"), c("A", "B"))))
  expect_identical(W[, c(3, 1)], cbind(Z, X))
  expect_identical(densify(W[c("a", "d", "b"), c("Z", "X")], 1, 1:2),
                   SimpleList(
                     list(list(Z = matrix(as.integer(c(1, 4, 2, 6, 9, 7)) + 1000L,
                                          ncol = 2,
                                          dimnames = list(c("a", "d", "b"),
                                                          c("A", "B"))),
                               X = matrix(as.integer(c(1, 4, 2, 6, 9, 7)),
                                          ncol = 2,
                                          dimnames = list(c("a", "d", "b"),
                                                          c("A", "B")))))))
})

test_that("[,SimplelistSparseAssays-method errors (out of bounds indices)", {
  msg <- "subscript contains NAs or out-of-bounds indices"
  expect_error(X[nrow(X) + 1, ], msg)
  expect_error(X["bad_idx", ], msg)
  expect_error(X[nrow(x) + 1, 1], msg)
  expect_error(X["bad_idx", 1], msg)
  expect_error(X[, ncol(x) + 1], msg)
  expect_error(X[, "bad_idx"], "subscript contains invalid names")
  expect_error(X[nrow(x) + 1, ncol(x) + 1], msg)
  expect_error(X["bad_idx", "another_bad_idx"], msg)
})

test_that("[,SimplelistSparseAssays-method warning (drop)", {
  expect_warning(X[1, , drop = TRUE],
                 "'drop' ignored '\\[,SimpleListSparseAssays,ANY,ANY-method'")
})

test_that("[<-,SimplelistSparseAssays-method errors (no index)", {
  # NOTE: While this might properly be considered a bug, subset replacing
  #       SparseAssays should be consistent with subsetting Assays in that at
  #       least one of i or j must be specified.
  expect_error(x[] <- x, "object 'fun' not found")
  expect_error(as(x, "ShallowSimpleListAssays")[] <-
                 as(x, "ShallowSimpleListAssays"))
})

test_that("[<-,SimplelistSparseAssays-method works (i only)", {
  # Test subsetting on i
  X <- x
  X[1:10, ] <- x[1:10, ]
  expect_that(x, not(is_identical_to(X)))
  expect_true(identical_SparseAssays(x, X))
  Y <- x
  Y[2, ] <- x[1, ]
  expect_identical(x[1, ], Y[2, ])
  expect_true(identical_SparseAssays(x[1, ], Y[2, ]))
})

test_that("[<-,SimplelistSparseAssays-method works (j only)", {
  # SparseAssays objects should be identical when only j is specified in subset
  # replacement
  X <- x
  X[, 1] <- x[, 1]
  expect_identical(x, X)
  Y <- x
  Y[, 2] <- x[, 1]
  # Sample names differ, so shouldn't be identical
  expect_that(x[, 1], not(is_identical_to(Y[, 2])))
  # Sample names differ, so shouldn't be identical
  expect_false(identical_SparseAssays(x[, 1], Y[, 2]))
  # But densified data should be identical
  expect_identical(densify(x[, 1], 1, 1)[[1L]][[1L]],
                   densify(Y[, 2], 1, 1)[[1L]][[1L]])
})

test_that("[<-,SimplelistSparseAssays-method works (i and j)", {
  X <- x
  X[1:7, 2] <- x[1:7, 2]
  expect_that(x, not(is_identical_to(X)))
  expect_true(identical_SparseAssays(x, X))
  Y <- x
  Y[2, 2] <- x[1, 1]
  # Sample names differ, so shouldn't be identical
  expect_that(x[1, 1], not(is_identical_to(Y[2, 2])))
  # Sample names differ, so shouldn't be identical
  expect_false(identical_SparseAssays(x[1, 1], Y[2, 2]))
  # But densified data should be identical
  expect_identical(densify(x[1, 1], 1, 1)[[1L]][[1L]],
                   densify(Y[2, 2], 1, 1)[[1L]][[1L]])
})

test_that("[<-,SimplelistSparseAssays-method works (non-contiguous i and j)", {
  X2 <- X
  x2 <- matrix(as.integer(c(11, 14, 12, 16, 19, 7)),
               ncol = 2, dimnames = list(c("a", "d", "b"), c("A", "B")))
  X2[c(1, 4, 2), ] <- SparseAssays(unname(x2))
  expect_identical(densify(X2[c(1, 4, 2), ], 1, 1)[[1L]][[1L]], x2)
  W2 <- W
  W2[, c(3, 1)] <- W[, c(1, 3)]
  # Sample names differ, so shouldn't be identical
  expect_that(W2, not(is_identical_to(W[, c(3, 1)])))
  # Sample names differ, so shouldn't be identical
  expect_false(identical_SparseAssays(W2, W[, c(3, 1)]))
  # But densified data should be identical
  for (i in seq_len(ncol(W))) {
    expect_identical(densify(W2, 1, 1:3)[[1L]][[i]],
                   densify(W[, 3:1], 1, 1:3)[[1L]][[i]])
  }
})

test_that("[<-,SimplelistSparseAssays-method works (non-numeric i and j)", {
  X2 <- X
  x2 <- matrix(as.integer(c(11, 14, 12, 16, 19, 7)),
               ncol = 2, dimnames = list(c("a", "d", "b"), c("A", "B")))
  X2[c("a", "d", "b"), ] <- SparseAssays(unname(x2))
  expect_identical(densify(X2[c("a", "d", "b"), ], 1, 1)[[1L]][[1L]],
                   x2)
  W2 <- W
  W2[, c("Z", "X")] <- W[, c("X", "Z")]
  # Sample names differ, so shouldn't be identical
  expect_that(W2, not(is_identical_to(W[, c("Z", "X")])))
  # Sample names differ, so shouldn't be identical
  expect_false(identical_SparseAssays(W2, W[, c("Z", "X")]))
  # But densified data should be identical
  for (i in seq_len(ncol(W))) {
    expect_identical(densify(W2, 1, 1:3)[[1L]][[i]],
                     densify(W[, 3:1], 1, 1:3)[[1L]][[i]])
  }
})

test_that("[<-,SimplelistSparseAssays-method errors (out of bounds indices)", {
  msg <- "subscript out of bounds"
  XX <- X
  expect_error(XX[nrow(X) + 1, ] <- X[1, ], msg)
  XX <- X
  expect_error(XX["bad_idx", ] <- X[1, ], msg)
  XX <- X
  expect_error(XX[nrow(X) + 1, 1] <- X[1, 1], msg)
  XX <- X
  expect_error(XX["bad_idx", 1] <- X[1, 1])
  XX <- X
  expect_error(XX[, ncol(X) + 1] <- X[, 1], msg)
  XX <- X
  expect_error(XX[, "bad_idx"] <- X[, 1], msg)
  XX <- X
  expect_error(XX[1, ncol(X) + 1] <- X[1, 1], msg)
  XX <- X
  expect_error(XX[1, "bad_idx"] <- X[1, 1], msg)
})

test_that("[<-,SimplelistSparseAssays-method errors (incorrectly dimensioned value)", {
  msg <- "number of items to replace is not a multiple of replacement length"
  expect_error(W[1, ] <- X[1, ], msg)
  expect_error(X[1, ] <- X[1:2, ], msg)
  expect_error(W[, 1] <- W[, 1:2], msg)
  expect_error(W[1, 1] <- W[1:2, 1:2], msg)
})

# TODO: Test rest of R/SimpleListSparseAssays-class.R
