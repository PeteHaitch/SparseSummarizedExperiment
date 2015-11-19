### =========================================================================
### Apply a function to a SparseAssays object in an [bp]lapply-ike manner
### -------------------------------------------------------------------------

setMethod("saapply", "SimpleListSparseAssays",
          function(X, FUN, densify = TRUE) {
            lapply(X, function(sparse_assay) {
              lapply(sparse_assay, function(sample) {
                if (densify) {
                  sample <- .densify.SimpleListSparseAssays.sample(sample)
                }
                FUN(sample)
              })
            }, FUN = FUN, densify = densify)
          }
)


# setMethod("saapply", "SimpleListSparseAssays",
#           function(X, FUN, densify = TRUE, ..., BPREDO = list(),
#                    BPPARAM = bpparam()) {
#             # TODO: How to capture dots?
#             dots <- list(...)
#             # TODO: How to pass down BPREDO and BPPARAM in a sensible fashion
#             bplapply(X, function(sparse_assay, FUN, densify, dots, BPREDO,
#                                  BPPARAM) {
#               bplapply(sparse_assay, function(sample, FUN, densify, dots,
#                                               BPREDO, BPPARAM) {
#                 if (densify) {
#                   sample <- .densify.SimpleListSparseAssays.sample(sample)
#                 }
#                 FUN(sample, ...)
#               }, FUN = FUN, densify = densify, dots = dots, BPREDO = BPREDO,
#               BPPARAM = BPPARAM)
#             }, FUN = FUN, densify = densify, dots = dots, BPREDO = BPREDO,
#             BPPARAM = BPPARAM)
#           }
# )

# TODO: Import BiocParallel
