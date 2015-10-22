#' An experimental implementation of a container modeled on the
#' SummarizedExperiment class for storing sparse genomic data.
#'
#' \pkg{SparseSummarizedExperiment} defines general purpose containers for
#' storing sparse genomic data. It aims to provide functionality for sparse
#' genomic data analogous to that available in the SummarizedExperiment
#' Bioconductor package.
#'
#' @docType package
#' @name SparseSummarizedExperiment-package
#'
# TODO: Try to be more selective about imports of methods and BiocGenerics and
#       include within the .R files
#' @importFrom methods setClass setGeneric setMethod callNextMethod validObject
# NOTE: Import the entire BiocGenerics package for simplicity
#' @import BiocGenerics
#'
# Below this line are old imports that I'm in the process of moving to the
# relevant .R files
# TODO: What's the correct way to import the C-routines from digest (if digest
#       is still required)?
NULL

# To avoid WARNINGs about "Undefined global functions or variables" in
# R CMD check
#' @importFrom utils globalVariables
globalVariables(c(".myI", ".myMap"))
