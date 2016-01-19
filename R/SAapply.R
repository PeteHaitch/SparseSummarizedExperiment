### =========================================================================
### Apply a function to a SparseAssays object in an bplapply-ike manner
### -------------------------------------------------------------------------

# TODO: How to pass down BPREDO and BPPARAM in a sensible fashion
# TODO: Only SerialParam backend is currently working; why?
# TODO: Should SAapply() return a list or a SimpleList when sparsify = FALSE?
#' @importFrom methods setMethod
#' @importFrom S4Vectors isTRUEorFALSE
#' @importFrom BiocParallel bplapply bpparam SerialParam
#'
#' @export
setMethod("SAapply", "SimpleListSparseAssays",
          # function(X, FUN, densify = TRUE, sparsify = !densify,
          #          withRownames = TRUE, ..., BPREDO = list(),
          #          BPPARAM = bpparam()) {
          function(X, FUN, densify = TRUE, sparsify = !densify,
                   withRownames = TRUE, ..., BPREDO = list(),
                   BPPARAM = SerialParam()) {
            if (!isTRUEorFALSE(densify) || !isTRUEorFALSE(sparsify)) {
              stop("'densify' and 'sparsify' must be TRUE or FALSE")
            }

            # NOTE: Can't pass 'densify', 'sparsify', and 'withRownames' via
            #       '...' in bplapply() because it is not an argument to 'FUN',
            #       so need to use a series of nested dummy functions.
            val <- bplapply(X, FUN = function(sparse_assay, densify, sparsify,
                                              withRownames, fun, ...) {
              bplapply(sparse_assay, function(sample, densify, sparsify,
                                              withRownames, fun, ...) {
                if (densify) {
                  sample <- .densify.SimpleListSparseAssays.sample(sample,
                                                                   withRownames)
                }
                val <- fun(sample, ...)
                if (!densify && sparsify) {
                  return(val)
                } else if (!densify && !sparsify) {
                  return(.densify.SimpleListSparseAssays.sample(val,
                                                                withRownames))
                } else if (densify && !sparsify) {
                  return(val)
                } else if (densify && sparsify) {
                  return(sparsify(val, "SimpleList"))
                }
              }, densify = densify, sparsify = sparsify,
              withRownames = withRownames, fun = FUN, ..., BPREDO = BPREDO,
              BPPARAM = BPPARAM)
            }, densify = densify, sparsify = sparsify,
            withRownames = withRownames, fun = FUN, ..., BPREDO = BPREDO,
            BPPARAM = BPPARAM)

            if (sparsify) {
              return(SparseAssays(val))
            }
            val
          }
)
