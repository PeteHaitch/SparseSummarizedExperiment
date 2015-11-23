context("Functionality that probably belongs in the SummarizedExperiment
        package")
test_that("combine,SummarizedExperiment0,SummarizedExperiment0-method works", {
  se0 <- se0[1:100, ]
  expect_identical(combine(se0, SummarizedExperiment::SummarizedExperiment()),
                   se0)
  expect_identical(combine(SummarizedExperiment::SummarizedExperiment(),
                           se0), se0)
  # Can't expect_identical() on SummarizedExperiment0 objects because assays
  # slot is a reference class.
  expect_equal(combine(se0[1:40, 1:4], se0[30:100, 3:6]), se0)
  expect_equal(combine(se0[1:40, ], se0[90:100, ]), se0[c(1:40, 90:100)])
  expect_equal(combine(se0[, 2], se0[, 2:4]), se0[, 2:4])
  se0_unnamed <- se0
  names(se0_unnamed) <- NULL
  expect_error(combine(se0_unnamed[, 1], se0_unnamed[, 2]),
               paste0("Cannot combine 'SummarizedExperiment0' objects with ",
                      "NULL 'names'"))
  se0_dupnames <- se0
  names(se0_dupnames) <- rep("A", nrow(se0_dupnames))
  expect_error(combine(se0_dupnames[, 2], se0_dupnames[, 1]),
               paste0("Cannot combine 'SummarizedExperiment0' objects with ",
                      "internal duplicate 'names'"))
  se0_nullcn <- se0
  colnames(se0_nullcn) <- NULL
  expect_error(combine(se0_nullcn[, 1], se0_nullcn[, 2]),
               paste0("Cannot combine 'SummarizedExperiment0' objects with ",
                      "NULL 'colnames'"))
  expect_error(combine(se0, rse),
               paste0("Cannot combine 'SummarizedExperiment0' and ",
                      "'RangedSummarizedExperiment' objects because only one ",
                      "of these has 'rowRanges'"))
})

test_that("combine,RangedSummarizedExperiment,RangedSummarizedExperiment0-method works", {
  rse <- rse[1:100, ]
  expect_identical(combine(rse, SummarizedExperiment::SummarizedExperiment()),
                   rse)
  expect_identical(combine(SummarizedExperiment::SummarizedExperiment(),
                           rse), rse)
  # Can't expect_identical() on SummarizedExperiment0 objects because assays
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


  expect_error(combine(rse, se0),
               paste0("Cannot combine 'RangedSummarizedExperiment' and ",
                      "'SummarizedExperiment0' objects because only one of ",
                      "these has 'rowRanges'"))
})
