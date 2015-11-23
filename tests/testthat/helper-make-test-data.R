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
