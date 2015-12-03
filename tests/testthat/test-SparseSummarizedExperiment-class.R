context("SparseSummarizedExperiment class")

test_that("SparseSummarizedExperiment class definition is correct", {
  expect_true(extends("SparseSummarizedExperiment", "SummarizedExperiment0"))
  expect_identical(slotNames("SparseSummarizedExperiment"),
                   c("sparseAssays", "colData", "assays", "NAMES",
                     "elementMetadata", "metadata"))
})

test_that(".valid.SSE.sparseAssays_class() and above work", {
  expect_null(.valid.SSE.sparseAssays_class(sse))
  expect_null(.valid.SSE(sse))
  expect_true(validObject(sse))
})

test_that(".valid.SSE.sparseAssays_nrow() and above", {
  v <- sse
  v@sparseAssays <- SparseAssays()
  expect_null(.valid.SSE.sparseAssays_nrow(v))
  expect_null(.valid.SSE.sparseAssays_dim(v))
  expect_null(.valid.SSE(v))
  expect_true(validObject(v))
  v <- sse
  v@sparseAssays <- sparseAssays(v)[1, ]
  msg <- "'sparseAssays' nrow differs from 'mcols' nrow"
  expect_identical(.valid.SSE.sparseAssays_nrow(v), msg)
  expect_identical(.valid.SSE.sparseAssays_dim(v), msg)
  expect_identical(.valid.SSE(v), msg)
  expect_error(validObject(v), msg)
})

test_that(".valid.SSE.sparseAssays_ncol() and above work", {
  v <- sse
  v@sparseAssays <- SparseAssays()
  expect_null(.valid.SSE.sparseAssays_ncol(v))
  expect_null(.valid.SSE.sparseAssays_dim(v))
  expect_null(.valid.SSE(v))
  expect_true(validObject(v))
  v <- sse
  v@sparseAssays <- sparseAssays(v)[, 1]
  msg <- "'sparseAssays' ncol differs from 'colData' nrow"
  expect_identical(.valid.SSE.sparseAssays_ncol(v), msg)
  expect_identical(.valid.SSE.sparseAssays_dim(v), msg)
  expect_identical(.valid.SSE(v), msg)
  expect_error(validObject(v), msg)
})

test_that("makeSEFromSSE works", {
  se <- makeSEFromSSE(sse)
  expect_true(SSE_identical_to_SE(sse, se))
  # TODO: RSE to RSSE coercion (should go in RSSE tests file)
})

test_that("Implicit coercion from SSE to SE0 works", {
  se0 <- as(sse, "SummarizedExperiment0")
  # Check all slots except sparseAssays, which is dropped in the coercion
  slot_names <- slotNames(sse)
  slot_names <- grep("sparseAssays", slot_names, value = TRUE, invert = TRUE)
  expect_true(all(sapply(slot_names, function(sn, x, y) {
    identical(slot(x, sn), slot(y, sn))
  }, x = sse, y = se0)))
})

test_that("Implicit coercion from SSE to RSE works", {
  rse <- as(sse, "RangedSummarizedExperiment")
  # Check all slots except sparseAssays, which is dropped in the coercion,
  # and rowRanges and NAMES, which are check separately since the SSE object
  # has no rowRanges and NAMES are transferred to the rowRanges slot in the
  # RSE object, are identical in the SSE and RSE.
  slot_names <- slotNames(sse)
  slot_names <- grep("sparseAssays|rowRanges|NAMES", slot_names, value = TRUE,
                     invert = TRUE)
  expect_true(all(sapply(slot_names, function(sn, x, y) {
    identical(slot(x, sn), slot(y, sn))
  }, x = sse, y = rse)))
  rr <- slot(rse, "rowRanges")
  expect_identical(length(rr), length(sse))
  expect_identical(names(rr), names(sse))
  expect_identical(names(rr), slot(sse, "NAMES"))
})

test_that("'Implicit' coercion from SSE to RSSE works", {
  rsse <- as(sse, "RangedSparseSummarizedExperiment")
  # Check all slots except rowRanges and NAMES, which are check separately
  # since the SSE object has no rowRanges and NAMES are transferred to the
  # rowRanges slot in the RSE object, are identical in the SSE and RSSE.
  slot_names <- slotNames(sse)
  slot_names <- grep("rowRanges|NAMES", slot_names, value = TRUE,
                     invert = TRUE)
  expect_true(all(sapply(slot_names, function(sn, x, y) {
    identical(slot(x, sn), slot(y, sn))
  }, x = sse, y = rsse)))
  rr <- slot(rsse, "rowRanges")
  expect_identical(length(rr), length(rsse))
  expect_identical(names(rr), names(rsse))
})

test_that("sparseAssays,SparseSummarizedExperiment-method works", {

})

test_that("sparseAssays<-,SparseSummarizedExperiment,SparseAssays-method works", {

})

test_that("sparseAssays<-,SparseSummarizedExperiment,SimpleList-method works", {

})

test_that("sparseAssays<-,SparseSummarizedExperiment,list-method works", {

})

test_that("sparseAssay,SparseSummarizedExperiment,missing-method", {

})

test_that("sparseAssay,SparseSummarizedExperiment,numeric-method", {

})

test_that("sparseAssay,SparseSummarizedExperiment,character-method", {

})

test_that("sparseAssay<-,SparseSummarizedExperiment,missing,SimpleList-method", {

})

test_that("sparseAssay<-,SparseSummarizedExperiment,numeric,SimpleList-method", {

})

test_that("sparseAssay<-,SparseSummarizedExperiment,character,SimpleList-method", {

})

test_that("sparseAssayNames,SparseSummarizedExperiment-method works", {

})

test_that("sparseAssayNames<-,SparseSummarizedExperiment,character-method works", {

})

test_that("[,SparseSummarizedExperiment-method works", {

})

test_that("[<-,SparseSummarizedExperiment-method works", {

})

test_that("rbind,SparseSummarizedExperiment-method works", {

})

test_that("cbind,SparseSummarizedExperiment-method works", {

})

test_that("combine,SparseSummarizedExperiment-method works", {

})
