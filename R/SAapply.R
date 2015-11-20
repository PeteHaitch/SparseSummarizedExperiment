### =========================================================================
### Apply a function to a SparseAssays object in an bplapply-ike manner
### -------------------------------------------------------------------------

# TODO: How to pass down BPREDO and BPPARAM in a sensible fashion
# TODO: Only SerialParam backend is currently working; why?
#' @importFrom methods setMethod
#' @importFrom BiocParallel bplapply bpparam SerialParam
#'
#' @export
setMethod("SAapply", "SimpleListSparseAssays",
          # function(X, FUN, densify = TRUE, ..., BPREDO = list(),
          #          BPPARAM = bpparam()) {
          function(X, FUN, densify = TRUE, ..., BPREDO = list(),
                   BPPARAM = SerialParam()) {

            # NOTE: Can't pass 'densify' via '...' in bplapply() because it is
            #       not an argument to 'FUN', so need to use a series of nested
            #       dummy functions.
            bplapply(X, FUN = function(sparse_assay, densify, fun, ...) {
              bplapply(sparse_assay, function(sample, densify, fun, ...) {
                if (densify) {
                  sample <- .densify.SimpleListSparseAssays.sample(sample)
                }
                fun(sample, ...)
              }, densify = densify, fun = FUN, ..., BPREDO = BPREDO,
              BPPARAM = BPPARAM)
            }, densify = densify, fun = FUN, ..., BPREDO = BPREDO,
            BPPARAM = BPPARAM)
          }
)
