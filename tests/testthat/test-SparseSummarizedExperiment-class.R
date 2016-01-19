context("SparseSummarizedExperiment class")

test_that("SparseSummarizedExperiment class definition is correct", {
  expect_true(extends("SparseSummarizedExperiment", "SummarizedExperiment"))
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
  msg <- paste0("\n  nb of rows in 'sparseAssays' (1) must equal nb of rows in ",
                "'rowData' (10000)")
  expect_identical(.valid.SSE.sparseAssays_nrow(v), msg)
  expect_identical(.valid.SSE.sparseAssays_dim(v), msg)
  expect_identical(.valid.SSE(v), msg)
  msg_escaped <- gsub("\n  ", "",
                      gsub("\\)", "\\\\)", gsub("\\(", "\\\\(", msg)))
  expect_error(validObject(v), msg_escaped)
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

test_that("Implicit coercion from SSE to SE works", {
  se <- as(sse, "SummarizedExperiment")
  # Check all slots except sparseAssays, which is dropped in the coercion
  slot_names <- slotNames(sse)
  slot_names <- grep("sparseAssays", slot_names, value = TRUE, invert = TRUE)
  expect_true(all(sapply(slot_names, function(sn, x, y) {
    identical(slot(x, sn), slot(y, sn))
  }, x = sse, y = se)))
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

test_that("sparseAssays,SparseSummarizedExperiment-method is compatible with generic", {
  generic <- getGeneric("sparseAssays")
  method <- getMethod("sparseAssays", "SparseSummarizedExperiment")
  expect_identical(generic@signature, "x")
  expect_identical(formals(generic@.Data), formals(method@.Data))
})

test_that("sparseAssays,SparseSummarizedExperiment-method works", {
  expect_true(
    identical_SparseAssays(sparseAssays(sse, withDimnames = FALSE),
                           sse@sparseAssays))
  expect_identical(
    .dimnames.SimpleListSparseAssays(sparseAssays(sse, withDimnames = TRUE)),
    dimnames(sse))
})

test_that("sparseAssays<-,SparseSummarizedExperiment,SparseAssays-method is compatible with generic", {
  generic <- getGeneric("sparseAssays<-")
  method <- getMethod("sparseAssays<-",
                      c("SparseSummarizedExperiment", "SparseAssays"))
  expect_identical(generic@signature, c("x", "value"))
  expect_identical(formals(generic@.Data), formals(method@.Data))
})

test_that("sparseAssays<-,SparseSummarizedExperiment,SparseAssays-method works", {
  SSE <- SparseSummarizedExperiment(x)
  sa1 <- SparseAssays(SimpleList(sa1 = as(x, "SimpleList")[[1]]))
  SSE2 <- SparseSummarizedExperiment(sa1)
  SSE1 <- SSE
  sparseAssays(SSE1) <- sa1
  expect_true(identical_SparseSummarizedExperiment(SSE1, SSE2))
})

test_that("sparseAssays<-,SparseSummarizedExperiment,SimpleList-method is compatible with generic", {
  generic <- getGeneric("sparseAssays<-")
  method <- getMethod("sparseAssays<-",
                      c("SparseSummarizedExperiment", "SimpleList"))
  expect_identical(generic@signature, c("x", "value"))
  expect_identical(formals(generic@.Data), formals(method@.Data))
})

test_that("sparseAssays<-,SparseSummarizedExperiment,SimpleList-method works", {
  SSE <- SparseSummarizedExperiment(x)
  sa1 <- SparseAssays(SimpleList(sa1 = as(x, "SimpleList")[[1]]))
  SSE2 <- SparseSummarizedExperiment(sa1)
  SSE1 <- SSE
  sparseAssays(SSE1) <- as(sa1, "SimpleList")
  expect_true(identical_SparseSummarizedExperiment(SSE1, SSE2))
})

test_that("sparseAssays<-,SparseSummarizedExperiment,list-method is compatible with generic", {
  generic <- getGeneric("sparseAssays<-")
  method <- getMethod("sparseAssays<-", c("SparseSummarizedExperiment", "list"))
  expect_identical(generic@signature, c("x", "value"))
  expect_identical(formals(generic@.Data), formals(method@.Data))
})

test_that("sparseAssays<-,SparseSummarizedExperiment,list-method works", {
  SSE <- SparseSummarizedExperiment(x)
  sa1 <- SparseAssays(SimpleList(sa1 = as(x, "SimpleList")[[1]]))
  SSE2 <- SparseSummarizedExperiment(sa1)
  SSE1 <- SSE
  sparseAssays(SSE1) <- as(sa1, "list")
  expect_true(identical_SparseSummarizedExperiment(SSE1, SSE2))
})

test_that("sparseAssay,SparseSummarizedExperiment,missing-method is compatible with generic", {
  generic <- getGeneric("sparseAssay")
  method <- getMethod("sparseAssay",
                      c("SparseSummarizedExperiment", "missing"))
  expect_identical(generic@signature, c("x", "i"))
  expect_identical(formals(generic@.Data), formals(method@.Data))
})

test_that("sparseAssay,SparseSummarizedExperiment,missing-method", {
  SSE <- SparseSummarizedExperiment(x)
  expect_true(identical_SparseAssays(sparseAssay(SSE),
                                     SparseAssays(SimpleList(sa1 = x1))))
  SSE0 <- SparseSummarizedExperiment()
  msg <- paste0("'sparseAssay\\(<SparseSummarizedExperiment>, i=\"missing\"",
                ", ...\\)' length\\(sparseAssays\\(",
                "<SparseSummarizedExperiment>\\)\\) is 0\n")
  expect_error(sparseAssay(SSE0), msg)
})

test_that("sparseAssay,SparseSummarizedExperiment,numeric-method is compatible with generic", {
  generic <- getGeneric("sparseAssay")
  method <- getMethod("sparseAssay",
                      c("SparseSummarizedExperiment", "numeric"))
  expect_identical(generic@signature, c("x", "i"))
  expect_identical(formals(generic@.Data), formals(method@.Data))
})

test_that("sparseAssay,SparseSummarizedExperiment,numeric-method", {
  SSE <- SparseSummarizedExperiment(x)
  expect_true(identical_SparseAssays(sparseAssay(SSE, 1),
                                     SparseAssays(SimpleList(sa1 = x1))))
  expect_true(identical_SparseAssays(sparseAssay(SSE, 2),
                                     SparseAssays(SimpleList(sa2 = x2))))
  msg <- paste0("'sparseAssay\\(<SparseSummarizedExperiment>, i=\"numeric\"",
                ", ...\\)' invalid subscript 'i'\n")
  expect_error(sparseAssay(SSE, 3), msg)
})

test_that("sparseAssay,SparseSummarizedExperiment,character-method is compatible with generic", {
  generic <- getGeneric("sparseAssay")
  method <- getMethod("sparseAssay",
                      c("SparseSummarizedExperiment", "character"))
  expect_identical(generic@signature, c("x", "i"))
  expect_identical(formals(generic@.Data), formals(method@.Data))
})

test_that("sparseAssay,SparseSummarizedExperiment,character-method", {
  SSE <- SparseSummarizedExperiment(x)
  expect_true(identical_SparseAssays(sparseAssay(SSE, "sa1"),
                                     SparseAssays(SimpleList(sa1 = x1))))
  expect_true(identical_SparseAssays(sparseAssay(SSE, "sa2"),
                                     SparseAssays(SimpleList(sa2 = x2))))
  msg <- paste0("'sparseAssay\\(<SparseSummarizedExperiment>, i=\"character\"",
                ", ...\\)' invalid subscript 'i'\n")
  expect_error(sparseAssay(SSE, "sa3"), msg)
})

test_that("sparseAssay<-,SparseSummarizedExperiment,missing,SimpleList-method is compatible with generic", {
  generic <- getGeneric("sparseAssay<-")
  method <- getMethod("sparseAssay<-",
                      c("SparseSummarizedExperiment", "missing", "SimpleList"))
  expect_identical(generic@signature, c("x", "i", "value"))
  expect_identical(formals(generic@.Data), formals(method@.Data))
})

test_that("sparseAssay<-,SparseSummarizedExperiment,missing,SimpleList-method", {
  SSE1 <- SparseSummarizedExperiment(X)
  X2 <- SparseAssays(matrix(100, ncol = ncol(X), nrow = nrow(Y),
                            dimnames = dimnames(X)))
  names(X2[[1L]]) <- names(X[[1L]])
  SSE2 <- SparseSummarizedExperiment(X2)
  sparseAssay(SSE1) <- as(X2, "SimpleList")[[1L]]
  expect_true(identical_SparseSummarizedExperiment(SSE1, SSE2))
  SSE0 <- SparseSummarizedExperiment()
  msg <- paste0("'sparseAssay\\(<SparseSummarizedExperiment>\\) <- value' ",
                "length\\(sparseAssays\\(",
                "<SparseSummarizedExperiment>\\)\\) is 0\n")
  expect_error(sparseAssay(SSE0) <- sparseAssays(SSE1), msg)
})

test_that("sparseAssay<-,SparseSummarizedExperiment,numeric,SimpleList-method is compatible with generic", {
  generic <- getGeneric("sparseAssay<-")
  method <- getMethod("sparseAssay<-",
                      c("SparseSummarizedExperiment", "numeric", "SimpleList"))
  expect_identical(generic@signature, c("x", "i", "value"))
  expect_identical(formals(generic@.Data), formals(method@.Data))
})

test_that("sparseAssay<-,SparseSummarizedExperiment,numeric,SimpleList-method", {
  SSE1 <- SparseSummarizedExperiment(X)
  X2 <- SparseAssays(matrix(100, ncol = ncol(X), nrow = nrow(Y),
                            dimnames = dimnames(X)))
  names(X2[[1L]]) <- names(X[[1L]])
  SSE2 <- SparseSummarizedExperiment(X2)
  sparseAssay(SSE1, 1) <- as(X2, "SimpleList")[[1L]]
  expect_true(identical_SparseSummarizedExperiment(SSE1, SSE2))
  # NOTE: Can add assay to end using sparseAssay<- method by effectively
  #       replacing a previously non-existant element.
  sparseAssay(SSE1, 2) <- as(X2, "SimpleList")[[1L]]
  expect_true(identical_SparseAssays(sparseAssay(SSE1, 1),
                                     sparseAssay(SSE1, 2)))
  # NOTE: But can't add assay beyond end using sparseAssay<- method (although
  #       the error message isn't very informative).
  msg <- paste0("All sparse assays of a 'SimpleListSparseAssays' object must ",
                "have an identical number of samples.")
  expect_error(sparseAssay(SSE1, 20) <- as(X2, "SimpleList")[[1L]], msg)
  # NOTE: names() must match
  msg <- paste0("colnames mismatch: all names\\(\\) of elements of a ",
                "'SimpleListSparseAssays' object must be identical.")
  expect_error(sparseAssay(SSE1, 1) <- as(Y, "SimpleList")[[1L]], msg)
})

test_that("sparseAssay<-,SparseSummarizedExperiment,character,SimpleList-method is compatible with generic", {
  generic <- getGeneric("sparseAssay<-")
  method <- getMethod("sparseAssay<-",
                      c("SparseSummarizedExperiment", "character", "SimpleList"))
  expect_identical(generic@signature, c("x", "i", "value"))
  expect_identical(formals(generic@.Data), formals(method@.Data))
})

test_that("sparseAssay<-,SparseSummarizedExperiment,character,SimpleList-method", {
  SSE1 <- SparseSummarizedExperiment(X)
  sparseAssayNames(SSE1) <- "kraken"
  X2 <- SparseAssays(matrix(100, ncol = ncol(X), nrow = nrow(Y),
                            dimnames = dimnames(X)))
  names(X2[[1L]]) <- names(X[[1L]])
  SSE2 <- SparseSummarizedExperiment(X2)
  sparseAssayNames(SSE2) <- "kraken"
  sparseAssay(SSE1, "kraken") <- as(X2, "SimpleList")[[1L]]
  expect_true(identical_SparseSummarizedExperiment(SSE1, SSE2))
  # NOTE: Can add sparse assay by name using sparseAssay<- method by
  #       effectively replacing a previously non-existant element.
  sparseAssay(SSE1, "godzilla") <- as(X2, "SimpleList")[[1L]]
  expect_true(identical_SparseAssays(unname(sparseAssay(SSE1, "kraken")),
                                     unname(sparseAssay(SSE1, "godzilla"))))
  # NOTE: names() must match
  msg <- paste0("colnames mismatch: all names\\(\\) of elements of a ",
                "'SimpleListSparseAssays' object must be identical.")
  expect_error(sparseAssay(SSE1, "kraken") <- as(Y, "SimpleList")[[1L]], msg)
})

test_that("sparseAssayNames,SparseSummarizedExperiment-method is compatible with generic", {
  generic <- getGeneric("sparseAssayNames")
  method <- getMethod("sparseAssayNames", "SparseSummarizedExperiment")
  expect_identical(generic@signature, "x")
  expect_identical(formals(generic@.Data), formals(method@.Data))
})

test_that("sparseAssayNames,SparseSummarizedExperiment-method works", {
  expect_identical(sparseAssayNames(sse), "sa1")
  expect_null(sparseAssayNames(SparseSummarizedExperiment()))
  expect_identical(sparseAssayNames(rsse), "sa1")
  expect_null(sparseAssayNames(
    SparseSummarizedExperiment(rowRanges = simGR(10L))))
})

test_that("sparseAssayNames<-,SparseSummarizedExperiment,character-method is compatible with generic", {
  generic <- getGeneric("sparseAssayNames<-")
  method <- getMethod("sparseAssayNames<-", c
                      ("SparseSummarizedExperiment", "character"))
  expect_identical(generic@signature, c("x", "value"))
  expect_identical(formals(generic@.Data), formals(method@.Data))
})

test_that("sparseAssayNames<-,SparseSummarizedExperiment,character-method works", {
  sparseAssayNames(sse) <- "kraken"
  expect_identical(sparseAssayNames(sse), "kraken")
  sparseAssayNames(rsse) <- "kraken"
  expect_identical(sparseAssayNames(rsse), "kraken")
  msg <- "invalid \\(NULL\\) left side of assignment"
  expect_error(sparseAssayNames(SparseSummarizedExperiment() <- "kraken"), msg)
  msg <- "'names' attribute \\[2\\] must be the same length as the vector"
  expect_error(sparseAssayNames(sse) <- c("kraken", "godzilla"), msg)
})

test_that("[,SparseSummarizedExperiment-method is compatible with generic", {
  generic <- getGeneric("[")
  method <- getMethod("[", "SparseSummarizedExperiment")
  expect_identical(generic@signature, c("x", "i", "j", "drop"))
  # NOTE: drop = TRUE for generic but drop = FALSE for method
  # expect_identical(formals(generic@.Data), formals(method@.Data))
})

test_that("[,SparseSummarizedExperiment-method works", {
  # UP TO HERE
})

test_that("[<-,SparseSummarizedExperiment-method is compatible with generic", {
  generic <- getGeneric("[<-")
  method <- getMethod("[<-",
                      c("SparseSummarizedExperiment", "ANY", "ANY",
                        "SparseSummarizedExperiment"))
  expect_identical(generic@signature, c("x", "i", "j", "value"))
  expect_identical(formals(generic@.Data), formals(method@.Data))
})

test_that("[<-,SparseSummarizedExperiment-method works", {

})

test_that("rbind,SparseSummarizedExperiment-method is compatible with generic", {
  generic <- getGeneric("rbind")
  method <- getMethod("rbind", "SparseSummarizedExperiment")
  expect_identical(generic@signature, "...")
  expect_identical(formals(generic@.Data), formals(method@.Data))
})

test_that("rbind,SparseSummarizedExperiment-method works", {

})

test_that("cbind,SparseSummarizedExperiment-method is compatible with generic", {
  generic <- getGeneric("cbind")
  method <- getMethod("cbind", "SparseSummarizedExperiment")
  expect_identical(generic@signature, "...")
  expect_identical(formals(generic@.Data), formals(method@.Data))
})

test_that("cbind,SparseSummarizedExperiment-method works", {

})

test_that("combine,SparseSummarizedExperiment-method is compatible with generic", {
  generic <- getGeneric("combine")
  method <- getMethod("combine",
                      c("SparseSummarizedExperiment",
                        "SparseSummarizedExperiment"))
  expect_identical(generic@signature, c("x", "y"))
  expect_identical(formals(generic@.Data), formals(method@.Data))
})

test_that("combine,SparseSummarizedExperiment-method works", {

})
