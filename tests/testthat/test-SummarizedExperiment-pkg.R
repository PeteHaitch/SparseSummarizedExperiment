context("Functionality that probably belongs in the SummarizedExperiment
        package")

# NOTE: This also coverage RangedSummarizedExperiment objects, since the method
#       is defined for these via inheritance to SummarizedExperiment.
test_that("combine,SummarizedExperiment,SummarizedExperiment-method is compatible with generic", {
  generic <- getGeneric("combine")
  method <- getMethod("combine",
                      c("SummarizedExperiment", "SummarizedExperiment"))
  expect_identical(generic@signature, c("x", "y"))
  expect_identical(formals(generic@.Data), formals(method@.Data))
})

test_that("combine,SummarizedExperiment,SummarizedExperiment-method works", {
  se <- se[1:100, ]
  expect_identical(combine(se, SummarizedExperiment::SummarizedExperiment()),
                   se)
  expect_identical(combine(SummarizedExperiment::SummarizedExperiment(),
                           se), se)
  # Can't expect_identical() on SummarizedExperiment objects because assays
  # slot is a reference class.
  expect_equal(combine(se[1:40, 1:4], se[30:100, 3:6]), se)
  expect_equal(combine(se[1:40, ], se[90:100, ]), se[c(1:40, 90:100)])
  expect_equal(combine(se[, 2], se[, 2:4]), se[, 2:4])
  se_unnamed <- se
  names(se_unnamed) <- NULL
  expect_error(combine(se_unnamed[, 1], se_unnamed[, 2]),
               paste0("Cannot combine 'SummarizedExperiment' objects with ",
                      "NULL 'names\\(\\)'"))
  se_dupnames <- se
  names(se_dupnames) <- rep("A", nrow(se_dupnames))
  expect_error(combine(se_dupnames[, 2], se_dupnames[, 1]),
               "'anyDuplicated\\(x\\)' must be 0 \\(FALSE\\)")
  se_nullcn <- se
  colnames(se_nullcn) <- NULL
  expect_error(combine(se_nullcn[, 1], se_nullcn[, 2]),
               paste0("Cannot combine 'SummarizedExperiment' objects with ",
                      "NULL 'colnames\\(\\)'"))
  expect_error(combine(se, rse),
               paste0("Cannot combine 'SummarizedExperiment' and ",
                      "'RangedSummarizedExperiment' objects because only one ",
                      "of these has a 'rowRanges' slot."))
})

test_that("combine,RangedSummarizedExperiment,RangedSummarizedExperiment-method works", {
  rse <- rse[1:100, ]
  expect_identical(combine(rse, SummarizedExperiment::SummarizedExperiment()),
                   rse)
  expect_identical(combine(SummarizedExperiment::SummarizedExperiment(),
                           rse), rse)
  # Can't expect_identical() on SummarizedExperiment objects because assays
  # slot is a reference class.
  expect_equal(combine(rse[1:40, 1:4], rse[30:100, 3:6]), rse)
  expect_equal(combine(rse[1:40, ], rse[90:100, ]), rse[c(1:40, 90:100)])
  expect_equal(combine(rse[, 2], rse[, 2:4]), rse[, 2:4])
  rse_dupranges <- rse
  rse_dupranges <- rbind(rse_dupranges[1:10, ], rse_dupranges[1:10, ])
  expect_error(combine(rse_dupranges[, 2], rse_dupranges[, 1]),
               "'any\\(duplicated\\(x\\)\\)' must be FALSE")
  rse_nullcn <- rse
  colnames(rse_nullcn) <- NULL
  expect_error(combine(rse_nullcn[, 1], rse_nullcn[, 2]),
               paste0("Cannot combine 'RangedSummarizedExperiment' objects ",
                      "with NULL 'colnames\\(\\)'"))
  rse_grl <- rse
  rowRanges(rse_grl) <- as(rowRanges(rse_grl), "GRangesList")
  expect_error(combine(rse_grl[1:40, 1:4], rse_grl[30:100, 3:6]),
               paste0("Cannot combine 'RangedSummarizedExperiment' objects ",
                      "with 'GRangesList'-based 'rowRanges'"))


  expect_error(combine(rse, se),
               paste0("Cannot combine 'RangedSummarizedExperiment' and ",
                      "'SummarizedExperiment' objects because only one of ",
                      "these has a 'rowRanges' slot."))
})

test_that("colnames are stripped from assays upon construction ", {
  # See https://stat.ethz.ch/pipermail/bioc-devel/2015-December/008410.html
  # This test is to keep track of the current behaviour in SummarizedExperiment
  m1 <- matrix(1:10, ncol = 2)
  m2 <- m1
  colnames(m2) <- c("A", "B")
  se1 <- SummarizedExperiment(m1, colData = DataFrame(row.names = c("A", "B")))
  se2 <- SummarizedExperiment(m2)
  se3 <- SummarizedExperiment(m2, colData = DataFrame(row.names = c("C", "D")))
  # colnames correctly set to c("A", "B") and stripped from assays
  expect_identical(colnames(se1), c("A", "B"))
  expect_null(colnames(se1@assays[[1L]]))
  # colnames correctly set to c("A", "B") set and but not stripped from assays
  expect_identical(colnames(se2), c("A", "B"))
  expect_identical(colnames(se2@assays[[1L]]), c("A", "B"))
  # colnames set to c("C", "D") (without warning about mismatch) and stripped
  # from assays
  expect_identical(colnames(se3), c("C", "D"))
  expect_null(colnames(se3@assays[[1L]]))
})
