### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Functions used to simulate data
###

#' Simulate data as a SimpleListSparseAssays object
#'
#' @param m Integer scalar giving the number of features.
#' @param n Integer scalar giving the number of samples.
#' @param d Numeric scalar in [0, 1] giving the proportion of non-missing data.
#' @param p A Integer scalar giving the ncol of the sparse assay.
#'
#' @return A SimpleListSparseAssays object
simSLSA <- function(m, n, d, p) {
  v <- replicate(n, {
    val <- matrix(NA_integer_, ncol = p, nrow = m)
    i <- sample(m, floor(d * m), replace = TRUE)
    val[i, ] <- rpois(floor(d * m) * p, lambda = 10)
    val
  }, simplify = FALSE)
  # Add sample names
  v <- mapply(function(vv, i) {
    tmp <- SparseAssays(vv)
    names(tmp[[1]]) <- paste0("s", i)
    tmp
  }, vv = v, i = seq_along(v))
  # Add sparse assay name
  v <- lapply(v, function(vv) {
    names(vv) <- "sa1"
    vv
  })
  # cbind (don't need to combine because know each element of v is a single
  # sample and has same dimensions by construction)
  do.call(cbind, v)
}

#' Simulate data as a GRanges object
#'
#' @param m Integer scalar giving the number of genomic ranges.
#'
#' @return A GRanges object.
simGR <- function(m) {
  GRanges(sample(paste0("chr", 1:22), m, replace = TRUE),
          IRanges(floor(runif(m, 1, 1e6)), width = 100),
          strand = sample(c("+", "-"), m, TRUE),
          feature_id = paste0("f", seq_len(m)))
}

#' Simulate data as a RangedSummarizedExperiment object
#'
#' @param m Integer scalar giving the number of genomic ranges.
#' @param n Integer scalar giving the number of samples (assumed/forced even).
#'
#' @return A RangedSummarizedExperiment object.
simRSE <- function(m, n) {
  # Require n to be even
  if (n %% 2 == 1) {
    n <- n + 1
  }
  counts <- matrix(floor(runif(m * n, 0, 1e4)), m)
  colData <- DataFrame(Treatment = rep(c("ChIP", "Input"), n / 2),
                       row.names = paste0("s", seq_len(n)))
  rowRanges <- simGR(m)
  SummarizedExperiment(assays = SimpleList(counts = counts),
                       rowRanges = rowRanges,
                       colData = colData)
}

#' Simulate data as a SummarizedExperiment0 object
#'
#' @param m Integer scalar giving the number of features.
#' @param n Integer scalar giving the number of samples (assumed/forced even).
#'
#' @return A SummarizedExperiment0 object.
simSE0 <- function(m, n) {
  as(simRSE(m, n), "SummarizedExperiment0")
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Helper functions
###

#' Test whether 2 SparseAssays are identical.
#'
#' Two SparseAssays are considered identical if their densified forms are
#' identical, i.e. identical_SparseAssays(x, y) is TRUE. Importantly, their
#' sparsified forms need be identical, i.e, identical(x, y) may be FALSE.
#'
#' @param x A SparseAssays object.
#' @param y A SparseAssays object.
#'
#' @return TRUE or FALSE
identical_SparseAssays <- function(x, y) {
  xx <- densify(x, seq_along(x), seq_len(ncol(x)))
  yy <- densify(y, seq_along(y), seq_len(ncol(y)))
  identical(xx, yy)
}

#' Test whether a SparseSummarizedExperiment is equivalent to a
#' SummarizedExperiment object.
#'
#' A SSE object is considered identical to a SE object if their colData,
#' NAMES/rowRanges, elementMetadata, metadata slots are identical, and the
#' SSE's sparseAssays + assays slot has equivalent data to the SE's assays slot.
#'
#' @param sse A SparseSummarizedExperiment or RangedSparseSummarizedExperiment
#'        object.
#' @param se A SummarizedExperiment0 or RangedSummarizedExperiment object.
#'
#' @return TRUE or FALSE
SSE_identical_to_SE <- function(sse, se) {
  if (!identical(sse@colData, se@colData)) {
    return(FALSE)
  }
  if ("rowRanges" %in% c(slotNames(sse), slotNames(se))) {
    if (!identical(sse@rowRanges, se@rowRanges)) {
      return(FALSE)
    }
  }
  if (!identical(sse@NAMES, se@NAMES)) {
    return(FALSE)
  }
  if (!identical(sse@elementMetadata, se@elementMetadata)) {
    return(FALSE)
  }
  if (!identical(sse@metadata, se@metadata)) {
    return(FALSE)
  }
  # NOTE: Need to remove sample names from sparseAssays before doing comparison
  #       since this is what makeSEFromSSE() does.
  unnamed_sa <- sse@sparseAssays
  unnamed_sa <- endoapply(unnamed_sa, function(x) {
    names(x) <- NULL
    x
  })
  # UP TO HERE: second has colnames (sample names)
  if (!identical(as(se@assays, "SimpleList")[sparseAssayNames(sse)],
                 as(as(unnamed_sa, "ShallowSimpleListAssays"),
                    "SimpleList"))) {
    return(FALSE)
  }
  return(TRUE)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Simulated objects used in tests
###

# NOTE: test data are generated using pseudorandom numbers, so need to set the
#       seed to ensure reproducibility.
set.seed(666)

m <- 10000
n <- 6
d <- 0.7
p <- 8

rse <- simRSE(m, n)
rsse <- as(rse, "RangedSparseSummarizedExperiment")
sa <- simSLSA(m, n, d, p)
sparseAssays(rsse) <- sa
sse <- as(rsse, "SparseSummarizedExperiment")
names(sse) <- paste0("F", seq_len(nrow(sse)))
se0 <- as(rse, "SummarizedExperiment0")
names(se0) <- paste0("F", seq_len(nrow(se0)))

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

XX <- matrix(1:10, ncol = 2, dimnames = list(letters[1:5], LETTERS[1:2]))
X <- SparseAssays(XX)
names(X[[1L]]) <- "X"
YY <- matrix(101:110, ncol = 2, dimnames = list(letters[1:5], LETTERS[1:2]))
Y <- SparseAssays(YY)
names(Y[[1L]]) <- "Y"
ZZ <- matrix(1001:1010, ncol = 2, dimnames = list(letters[1:5], LETTERS[1:2]))
Z <- SparseAssays(ZZ)
names(Z[[1L]]) <- "Z"
W <- cbind(X, Y, Z)
