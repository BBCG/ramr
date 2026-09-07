# Merges, filters and outputs all genomic regions of a given \`GRanges\` object

\`getUniverse\` returns a \`GRanges\` object with all the genomic
regions in a data set, that can be used for AMR enrichment analysis

## Usage

``` r
getUniverse(data.ranges, merge.window = 300, min.cpgs = 7, min.width = 1)
```

## Arguments

- data.ranges:

  A \`GRanges\` object with genomic locations and corresponding beta
  values included as metadata.

- merge.window:

  A single integer \>= 1. All \`data.ranges\` genomic locations within
  this distance will be merged (the default: 300).

- min.cpgs:

  A single integer \>= 1. All genomic regions containing less than
  \`min.cpgs\` genomic locations are filtered out (the default: 7).

- min.width:

  A single integer \>= 1 (the default). Only regions with the width of
  at least \`min.width\` are returned.

## Value

The output is a \`GRanges\` object that contain all the genomic regions
in \`data.ranges\` object (in other words, all potential AMRs).

## Details

In the provided data set \`getUniverse\` merges and outputs all the
genomic regions that satisfy filtering criteria, thus creating a
\`GRanges\` object to be used as a reference set of genomic regions for
AMR enrichment analysis.

## See also

[`getAMR`](getAMR.md) for identification of AMRs,
[`plotAMR`](plotAMR.md) for plotting AMRs,
[`simulateAMR`](simulateAMR.md) and [`simulateData`](simulateData.md)
for the generation of simulated test data sets, and \`ramr\` vignettes
for the description of usage and sample data.

## Examples

``` r
  data(ramr)
  universe <- getUniverse(ramr.data, min.cpgs=5, merge.window=1000)
# \donttest{
  # identify AMRs
  amrs <- getAMR(
    data.ranges=ramr.data, compute="beta+binom",
    combine.min.cpgs=5, combine.window=1000, combine.threshold=1e-3
  )
#> Preprocessing data 
#> [0.057s]
#> Fitting beta distribution 
#> [0.035s]
#> Creating genomic ranges 
#> [0.012s]

  # AMR enrichment analysis using LOLA
  library(LOLA)
  # download LOLA region databases from http://databio.org/regiondb
  hg19.extdb.file <- system.file("LOLAExt", "hg19", package="LOLA")
  if (file.exists(hg19.extdb.file)) {
    hg19.extdb  <- loadRegionDB(hg19.extdb.file)
    runLOLA(amrs, universe, hg19.extdb, cores=1, redefineUserSets=TRUE)
  }
# }
```
