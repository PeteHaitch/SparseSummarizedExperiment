# Combining SummarizedExperiment objects
Peter Hickey  
15 October 2015  



# Motivation

I often find myself with multiple `SE` objects (I'm using `SE` as a shorthand for the `SummarizedExperiment0` and `RangedSummarizedExeriment` classes), each with different samples but possibly non-overlapping features/ranges. Currently, it is difficult to combine these objects;  `rbind()` can only combine objects with the same samples but different features/ranges and `cbind()` can only combine objects with the same features/ranges but different samples. I think it would be useful to have a "combine" method for `SE` objects that handles the situation where each object has different samples but with possibly non-overlapping features/ranges.

# Methods

## `combineSE()`

This is a first pass at addressing this need.


```r
suppressPackageStartupMessages(library(SummarizedExperiment))
# Need to set seed because example("RangedSumamrizedExperiment") and 
# example("SummarizedExperiment") both call runif().
set.seed(666)
```

`combineSE()` is the workhorse function for combining `SE` objects. There are a few things I want to point out:

1. `colnames` (which I think of as sample names) must be unique across all objects.
2. All objects must have the same `assays`.
3. If the objects are `SummarizedExperiment0` objects, then the `NAMES` slot must be non-NULL. This is because matching of features across objects is done using the `NAMES`.
4. The value of the `nomatch` argument is used to fill in missing values, e.g, where a sample does not have a value for that particular feature/range in that assay.
5. Combining the `elementMetadata` is pretty rough. For `SummarizedExperiment0` objects, combining `elementMetadata` is pretty inefficient since it calls the `unique,DataFrame-method`. Therefore the default is to drop the metadata columns (`use.mcols = FALSE`). There is probably a better method for combining the `elementMetadata` of `SE` objects.


```r
combineSE <- function(x, y, ..., nomatch = NA, use.mcols = FALSE) {
            
            args <- unname(list(x, y, ...))
            
            # Check that each object has unique colnames
            colnames <- unlist(lapply(args, colnames))
            if (anyDuplicated(colnames)) {
              stop(paste0("Cannot combine ", class(args[[1]]), " objects with ", 
                          "duplicate 'colnames'"))
            }
            
            # Check that each object has the same assays
            an <- lapply(args, assayNames)
            if (any(sapply(an, function(x, y) any(is.na(match(x, y))), 
                           y = an[[1]]))) {
              stop(paste0("All ", class(args[[1]]), " objects must have ", 
                          "identical 'assayNames'"))
            }
            
            # Combine rowRanges or NAMES slot
            if (is(args[[1L]], "RangedSummarizedExperiment")) {
              # NOTE: rowRanges(x) includes elementMetadata(x) since 
              #       elementMetadata(x) = elementMetadata(rowRanges(x)).
              #       This is the reason for the use.mcols if-else block.
              if (use.mcols) {
                all_rowRanges <- do.call(c, lapply(args, rowRanges))
              } else {
                all_rowRanges <- do.call(c, lapply(args, function(x) {
                  rr <- rowRanges(x)
                  elementMetadata(rr) <- NULL
                  rr
                  }))
              }
              rowRanges <- unique(all_rowRanges)
              nr <- length(rowRanges)
            } else {
              if (any(vapply(args, function(x) is.null(x@NAMES), logical(1)))) {
                stop("Cannot combine ", class(args[[1]]), " objects with ", 
                     "'NAMES' set to NULL")
              }
              all_NAMES <- do.call(c, lapply(args, function(x) x@NAMES))
              NAMES <- unique(all_NAMES)
              nr <- length(NAMES)
            }
            
            # Combine colData
            colData <- do.call(rbind, lapply(args, colData))
            
            # Combine elementMetadata
            if (use.mcols) {
              if (is(args[[1L]], "RangedSummarizedExperiment")) {
                elementMetadata <- mcols(rowRanges)
              } else {
                # IDEA: Create DataFrame with all_NAMES/all_rowRanges in one 
                #       column and elementMetadata in others, unique-ify, and 
                #       check that the number of unique rows equals nr.
                # WARNING: This will be slow for large SummarizedExperiment0 
                #          objects
                elementMetadata <- unique(
                  cbind(DataFrame(all_NAMES), 
                        do.call(rbind, lapply(args, elementMetadata))))
                # Now drop the all_NAMES column (assumes it is the first column)
                elementMetadata <- elementMetadata[, -c(1L), drop = FALSE]
                # Sanity check
                if (nrow(elementMetadata) != nr) {
                  stop(paste0("'elementMetadata' must match across ", 
                              class(args[[1]]), " objects"))
                }
              }
            } else {
              elementMetadata <- DataFrame()
              elementMetadata@nrows <- nr
            }
            
            # Create assays of the correct dimension (fill with 'nomatch')
            # First, create the empty combined assay using the appropriate 
            # storage.mode (guessed from the storage.mode of the assay in the 
            # first sample).
            nomatch <- lapply(seq_along(an[[1]]), function(i) {
              storage.mode(nomatch) <- storage.mode(assay(args[[1]], 
                                                          i, 
                                                          withDimnames = FALSE))
              nomatch
            })
            assays <- lapply(nomatch, function(nm) {
              matrix(nm, nrow = nr, ncol = length(colnames))
            })
            names(assays) <- an[[1]]
            
            # NOTE: I suspect that there are faster and more efficient ways to 
            # combine the assays, perhaps at the C-level.
            if (is(args[[1L]], "RangedSummarizedExperiment")) {
              for (j in seq_along(args)) {
                ol <- findOverlaps(args[[j]], rowRanges, type = "equal")
                for (i in seq_along(assays)) {
                  assays[[i]][subjectHits(ol),
                              match(colnames(args[[j]]), colnames)] <- 
                    assay(args[[j]], i, withDimnames = FALSE)
                }
              }
            } else {
              for (j in seq_along(args)) {
                ol <- match(args[[j]]@NAMES, NAMES)
                for (i in seq_along(assays)) {
                  assays[[i]][ol, match(colnames(args[[j]]), colnames)] <- 
                    assay(args[[j]], i, withDimnames = FALSE)
                }
              }
            }
            assays <- Assays(assays)
            
            # Combine metadata
            metadata <- do.call(c, lapply(args, metadata))
            
            if (is(args[[1L]], "RangedSummarizedExperiment")) {
              # No need to replace elementMetadata slot since it is part of 
              # rowRanges.
              BiocGenerics:::replaceSlots(args[[1L]], 
                                          rowRanges = rowRanges,
                                          colData = colData, 
                                          assays = assays,
                                          metadata = metadata)
            } else {
              BiocGenerics:::replaceSlots(args[[1L]],
                                          NAMES = NAMES,
                                          colData = colData,
                                          assays = assays,
                                          metadata = metadata,
                                          elementMetadata = elementMetadata)
            }
          }
```

## The need for a new generic

Unfortunately, the `BiocGenerics::combine` generic will not work for with `combineSE()` (I think it's because of the use of the `nomatch` and `use.mcols` arguments in `combineSE()`). Also, `BiocGenerics::combine` works recursively, which unnecessarily slows things down since `combineSE()` can work on all arguments in the one function call.


```r
BiocGenerics::combine
#> nonstandardGenericFunction for "combine" defined from package "BiocGenerics"
#> 
#> function (x, y, ...) 
#> {
#>     if (length(list(...)) > 0L) {
#>         combine(x, do.call(combine, list(y, ...)))
#>     }
#>     else {
#>         standardGeneric("combine")
#>     }
#> }
#> <environment: 0x7faa42f448d0>
#> Methods may be defined for arguments: x, y
#> Use  showMethods("combine")  for currently available ones.
```

I therefore define a new generic, which for now I call `combine2`


```r
setGeneric("combine2", function(x, y, ...) {
  standardGeneric("combine2")
})
#> [1] "combine2"

# NOTE: Also defined for RangedSummarizedExperiment via inheritance
setMethod("combine2", c("SummarizedExperiment0", "SummarizedExperiment0"),
          function(x, y, ..., nomatch = NA, use.mcols = FALSE) {
            combineSE(x, y, ..., nomatch = nomatch, use.mcols = use.mcols)
          }
)
#> [1] "combine2"
```

## Examples

I use data from the `SummarizedExperiment::SummarizedExeriment` and `SummarizedExperimentRangedSummarizedExperiment` man pages.

### `SummarizedExperiment0` objects


```r
# Get some example data
example("SummarizedExperiment", echo = FALSE)

# NOTE: SummarizedExperiment0 objects must have non-NULL NAMES slot in order 
#       to combine
se0@NAMES <- as.character(1:200)

# Create data to combine
A <- se0[1:4, "A"]
B <- se0[3:6, "B"]
C <- se0[5:8, "C"]
D <- se0[7:10, "D"]
E <- se0[9:12, "E"]
F <- se0[11:14, "F"]

# Sanity check: identical to cbind when given compatible arguments
all.equal(cbind(se0[, 1], se[, 2], se[, 3]),
          combine2(se0[, 1], se0[, 2], se0[, 3], use.mcols = TRUE))
#> [1] TRUE

# The default
x <- combine2(A, B, C, D, E, F)
x
#> class: SummarizedExperiment0 
#> dim: 14 6 
#> metadata(0):
#> assays(1): counts
#> rownames(14): 1 2 ... 13 14
#> metadata column names(0):
#> colnames(6): A B ... E F
#> colData names(1): Treatment
assay(x)
#>           A        B        C        D        E        F
#> 1  9.647809       NA       NA       NA       NA       NA
#> 2  8.280480       NA       NA       NA       NA       NA
#> 3  9.881258 6.543259       NA       NA       NA       NA
#> 4  8.301061 8.931229       NA       NA       NA       NA
#> 5        NA 9.879669 8.980169       NA       NA       NA
#> 6        NA 9.438820 8.575940       NA       NA       NA
#> 7        NA       NA 9.808847 8.695873       NA       NA
#> 8        NA       NA 9.727737 8.913906       NA       NA
#> 9        NA       NA       NA 6.484487 8.890558       NA
#> 10       NA       NA       NA 7.971683 9.186419       NA
#> 11       NA       NA       NA       NA 9.681241 9.880142
#> 12       NA       NA       NA       NA 9.351869 8.603473
#> 13       NA       NA       NA       NA       NA 9.238812
#> 14       NA       NA       NA       NA       NA 7.739453

# Using mcols
y <- combine2(A, B, C, D, E, F, use.mcols = TRUE)
y
#> class: SummarizedExperiment0 
#> dim: 14 6 
#> metadata(0):
#> assays(1): counts
#> rownames(14): 1 2 ... 13 14
#> metadata column names(1): feature_id
#> colnames(6): A B ... E F
#> colData names(1): Treatment
mcols(y)
#> DataFrame with 14 rows and 1 column
#>      feature_id
#>     <character>
#> 1         ID001
#> 2         ID002
#> 3         ID003
#> 4         ID004
#> 5         ID005
#> ...         ...
#> 10        ID010
#> 11        ID011
#> 12        ID012
#> 13        ID013
#> 14        ID014

# Using -99 for nomatch
z <- combine2(A, B, C, D, E, F, nomatch = -99)
assay(z)
#>             A          B          C          D          E          F
#> 1    9.647809 -99.000000 -99.000000 -99.000000 -99.000000 -99.000000
#> 2    8.280480 -99.000000 -99.000000 -99.000000 -99.000000 -99.000000
#> 3    9.881258   6.543259 -99.000000 -99.000000 -99.000000 -99.000000
#> 4    8.301061   8.931229 -99.000000 -99.000000 -99.000000 -99.000000
#> 5  -99.000000   9.879669   8.980169 -99.000000 -99.000000 -99.000000
#> 6  -99.000000   9.438820   8.575940 -99.000000 -99.000000 -99.000000
#> 7  -99.000000 -99.000000   9.808847   8.695873 -99.000000 -99.000000
#> 8  -99.000000 -99.000000   9.727737   8.913906 -99.000000 -99.000000
#> 9  -99.000000 -99.000000 -99.000000   6.484487   8.890558 -99.000000
#> 10 -99.000000 -99.000000 -99.000000   7.971683   9.186419 -99.000000
#> 11 -99.000000 -99.000000 -99.000000 -99.000000   9.681241   9.880142
#> 12 -99.000000 -99.000000 -99.000000 -99.000000   9.351869   8.603473
#> 13 -99.000000 -99.000000 -99.000000 -99.000000 -99.000000   9.238812
#> 14 -99.000000 -99.000000 -99.000000 -99.000000 -99.000000   7.739453
```

The error message is pretty ugly and uninformative if the `elementMetadata` aren't compatible across `SE`.


```r
# An ugly error due to incompatible elementMetadata
a <- A
elementMetadata(a) <- DataFrame(J = seq_len(nrow(a)))
combine2(a, B, use.mcols = TRUE)
#> Error in unique(cbind(DataFrame(all_NAMES), do.call(rbind, lapply(args, : error in evaluating the argument 'x' in selecting a method for function 'unique': Error in .Method(..., deparse.level = deparse.level) : 
#>   column names for arg 2 do not match those of first arg
#> Calls: cbind ... <Anonymous> -> standardGeneric -> eval -> eval -> eval -> .Method
```

### `RangedSummarizedExperiment` objects


```r
# Get some example data
example("RangedSummarizedExperiment", echo = FALSE)

# Create 
A <- rse[1:4, "A"]
B <- rse[3:6, "B"]
C <- rse[5:8, "C"]
D <- rse[7:10, "D"]
E <- rse[9:12, "E"]
F <- rse[11:14, "F"]

# Sanity check: identical to cbind when given compatible arguments
all.equal(cbind(rse[, 1], rse[, 2], rse[, 3]),
          combine2(rse[, 1], rse[, 2], rse[, 3], use.mcols = TRUE))
#> [1] TRUE

# The default
x <- combine2(A, B, C, D, E, F)
x
#> class: RangedSummarizedExperiment 
#> dim: 4 6 
#> metadata(0):
#> assays(1): counts
#> rownames: NULL
#> rowRanges metadata column names(0):
#> colnames(6): A B ... E F
#> colData names(1): Treatment
assay(x)
#>              A        B        C        D        E        F
#>  [1,] 8.753012       NA       NA       NA       NA       NA
#>  [2,] 9.450510       NA       NA       NA       NA       NA
#>  [3,] 8.653947 8.817888       NA       NA       NA       NA
#>  [4,] 9.725219 8.300851       NA       NA       NA       NA
#>  [5,]       NA 9.707305 8.185367       NA       NA       NA
#>  [6,]       NA 9.880137 7.923677       NA       NA       NA
#>  [7,]       NA       NA 7.840996 9.796866       NA       NA
#>  [8,]       NA       NA 8.674516 9.528341       NA       NA
#>  [9,]       NA       NA       NA 9.453123 8.651413       NA
#> [10,]       NA       NA       NA 9.532267 9.147541       NA
#> [11,]       NA       NA       NA       NA 9.119392 9.279716
#> [12,]       NA       NA       NA       NA 9.080164 9.566621
#> [13,]       NA       NA       NA       NA       NA 9.839014
#> [14,]       NA       NA       NA       NA       NA 9.704396

# Using mcols
y <- combine2(A, B, C, D, E, F, use.mcols = TRUE)
y
#> class: RangedSummarizedExperiment 
#> dim: 4 6 
#> metadata(0):
#> assays(1): counts
#> rownames: NULL
#> rowRanges metadata column names(1): feature_id
#> colnames(6): A B ... E F
#> colData names(1): Treatment
mcols(y)
#> DataFrame with 14 rows and 1 column
#>      feature_id
#>     <character>
#> 1         ID001
#> 2         ID002
#> 3         ID003
#> 4         ID004
#> 5         ID005
#> ...         ...
#> 10        ID010
#> 11        ID011
#> 12        ID012
#> 13        ID013
#> 14        ID014

# Using -99 for nomatch
z <- combine2(A, B, C, D, E, F, nomatch = -99)
assay(z)
#>                A          B          C          D          E          F
#>  [1,]   8.753012 -99.000000 -99.000000 -99.000000 -99.000000 -99.000000
#>  [2,]   9.450510 -99.000000 -99.000000 -99.000000 -99.000000 -99.000000
#>  [3,]   8.653947   8.817888 -99.000000 -99.000000 -99.000000 -99.000000
#>  [4,]   9.725219   8.300851 -99.000000 -99.000000 -99.000000 -99.000000
#>  [5,] -99.000000   9.707305   8.185367 -99.000000 -99.000000 -99.000000
#>  [6,] -99.000000   9.880137   7.923677 -99.000000 -99.000000 -99.000000
#>  [7,] -99.000000 -99.000000   7.840996   9.796866 -99.000000 -99.000000
#>  [8,] -99.000000 -99.000000   8.674516   9.528341 -99.000000 -99.000000
#>  [9,] -99.000000 -99.000000 -99.000000   9.453123   8.651413 -99.000000
#> [10,] -99.000000 -99.000000 -99.000000   9.532267   9.147541 -99.000000
#> [11,] -99.000000 -99.000000 -99.000000 -99.000000   9.119392   9.279716
#> [12,] -99.000000 -99.000000 -99.000000 -99.000000   9.080164   9.566621
#> [13,] -99.000000 -99.000000 -99.000000 -99.000000 -99.000000   9.839014
#> [14,] -99.000000 -99.000000 -99.000000 -99.000000 -99.000000   9.704396
```

The error message is pretty ugly and uninformative if the `elementMetadata` aren't compatible across `SE`.


```r
# An ugly error due to incompatible elementMetadata
a <- A
elementMetadata(a) <- DataFrame(J = seq_len(nrow(a)))
combine2(a, B, use.mcols = TRUE)
#> Error in .Method(..., deparse.level = deparse.level): column names for arg 2 do not match those of first arg
```

# Summary

The `combine2()` method for `SE` objects address my initial aim of being able to combine multiple `SE` objects when they have different samples but possibly non-overlapping features/ranges.

## Futher work

- A better name than `combine2`. Can we use `combine` (I couldn't get it to be compatible with the `BiocGenerics::combine` generic)?
- Unit tests
- Documentation
- A more informative error message if the user tries to combine `RangedSummarizedExperiment` and `SummarizedExperiment0` objects together.
- A more informative error message if `elementMetadata` aren't compatible across `SE` objects.
- [__low priority__] Extending `combine2` to work with objects that contain data from the same samples and possibly non-overlapping features/ranges? It will require a check that "overlapping" measurements are identical and doing something appropriate if they aren't.
- A general method for "combining" `SE` objects, e.g., you don't need to know whether you should be doing a `rbind()`/`cbind()`/`combine2()`, the "combining" method dispatches to the appropriate sub-method automagically.

# Session info


```r
devtools::session_info()
#> Session info --------------------------------------------------------------
#>  setting  value                                             
#>  version  R Under development (unstable) (2015-10-13 r69511)
#>  system   x86_64, darwin13.4.0                              
#>  ui       X11                                               
#>  language (EN)                                              
#>  collate  en_AU.UTF-8                                       
#>  tz       Australia/Melbourne                               
#>  date     2015-10-15
#> Packages ------------------------------------------------------------------
#>  package              * version date      
#>  Biobase              * 2.31.0  2015-10-14
#>  BiocGenerics         * 0.17.0  2015-10-14
#>  devtools               1.9.1   2015-09-11
#>  digest                 0.6.8   2014-12-31
#>  evaluate               0.8     2015-09-18
#>  formatR                1.2.1   2015-09-18
#>  GenomeInfoDb         * 1.7.0   2015-10-14
#>  GenomicRanges        * 1.21.32 2015-10-14
#>  htmltools              0.2.6   2014-09-08
#>  IRanges              * 2.4.0   2015-10-14
#>  knitr                  1.11    2015-08-14
#>  magrittr               1.5     2014-11-22
#>  memoise                0.2.1   2014-04-22
#>  rmarkdown              0.8.1   2015-10-10
#>  S4Vectors            * 0.9.0   2015-10-14
#>  stringi                0.5-5   2015-06-29
#>  stringr                1.0.0   2015-04-30
#>  SummarizedExperiment * 1.1.0   2015-10-14
#>  XVector                0.11.0  2015-10-14
#>  yaml                   2.1.13  2014-06-12
#>  zlibbioc               1.17.0  2015-10-14
#>  source                                      
#>  Bioconductor                                
#>  Bioconductor                                
#>  CRAN (R 3.3.0)                              
#>  CRAN (R 3.3.0)                              
#>  CRAN (R 3.3.0)                              
#>  CRAN (R 3.3.0)                              
#>  Bioconductor                                
#>  Bioconductor                                
#>  CRAN (R 3.3.0)                              
#>  Github (Bioconductor-mirror/IRanges@b0b5fff)
#>  CRAN (R 3.3.0)                              
#>  CRAN (R 3.3.0)                              
#>  CRAN (R 3.3.0)                              
#>  CRAN (R 3.3.0)                              
#>  Bioconductor                                
#>  CRAN (R 3.3.0)                              
#>  CRAN (R 3.3.0)                              
#>  Bioconductor                                
#>  Bioconductor                                
#>  CRAN (R 3.3.0)                              
#>  Bioconductor
```
