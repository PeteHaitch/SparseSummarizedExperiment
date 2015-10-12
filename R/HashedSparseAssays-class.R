### =========================================================================
### HashedSparseAssays objects
### -------------------------------------------------------------------------
###
### This is just an idea at this stage. Within each sparse assay, each sample
### would have a vector of keys (similar to the the 'map' elements in the
### SparseAssays object) that is the key in a hash table. Rather than than
### having a 'data' element for each sample, the HashedSparseAssays class would
### use a single hash table 'data' element per sparse assay.
