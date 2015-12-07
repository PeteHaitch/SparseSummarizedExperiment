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
                      "NULL 'names'"))
  se_dupnames <- se
  names(se_dupnames) <- rep("A", nrow(se_dupnames))
  expect_error(combine(se_dupnames[, 2], se_dupnames[, 1]),
               paste0("Cannot combine 'SummarizedExperiment' objects with ",
                      "internal duplicate 'names'"))
  se_nullcn <- se
  colnames(se_nullcn) <- NULL
  expect_error(combine(se_nullcn[, 1], se_nullcn[, 2]),
               paste0("Cannot combine 'SummarizedExperiment' objects with ",
                      "NULL 'colnames'"))
  expect_error(combine(se, rse),
               paste0("Cannot combine 'SummarizedExperiment' and ",
                      "'RangedSummarizedExperiment' objects because only one ",
                      "of these has 'rowRanges'"))
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
               paste0("Cannot combine 'RangedSummarizedExperiment' objects ",
                      "with internal duplicate 'rowRanges'"))
  rse_nullcn <- rse
  colnames(rse_nullcn) <- NULL
  expect_error(combine(rse_nullcn[, 1], rse_nullcn[, 2]),
               paste0("Cannot combine 'RangedSummarizedExperiment' objects ",
                      "with NULL 'colnames'"))
  rse_grl <- rse
  rowRanges(rse_grl) <- as(rowRanges(rse_grl), "GRangesList")
  expect_error(combine(rse_grl[1:40, 1:4], rse_grl[30:100, 3:6]),
               paste0("Cannot combine 'RangedSummarizedExperiment' objects ",
                      "with 'GRangesList'-based 'rowRanges'"))


  expect_error(combine(rse, se),
               paste0("Cannot combine 'RangedSummarizedExperiment' and ",
                      "'SummarizedExperiment' objects because only one of ",
                      "these has 'rowRanges'"))
})
