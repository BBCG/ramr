# \[OBSOLETE\] Search for aberrantly methylated regions

This function is fully functional but obsolete. It will remain a part of
the package for consistency, as it was used in \`ramr\` publication
([doi:10.1093/bioinformatics/btab586](https://doi.org/10.1093/bioinformatics/btab586)
). Please use faster and more capable [`getAMR`](getAMR.md) instead.

\`getAMR.obsolete\` returns a \`GRanges\` object with all the aberrantly
methylated regions (AMRs) for all samples in a data set.

## Usage

``` r
getAMR.obsolete(
  data.ranges,
  data.samples = NULL,
  ramr.method = c("IQR", "beta", "wbeta", "beinf"),
  iqr.cutoff = 5,
  pval.cutoff = 0.05,
  qval.cutoff = NULL,
  merge.window = 300,
  min.cpgs = 7,
  min.width = 1,
  exclude.range = NULL,
  cores = max(1, parallel::detectCores() - 1),
  verbose = TRUE,
  ...
)
```

## Arguments

- data.ranges:

  A \`GRanges\` object with genomic locations and corresponding beta
  values included as metadata.

- data.samples:

  A character vector with sample names (a subset of metadata column
  names). If \`NULL\` (the default), then all samples (metadata columns)
  are included in the analysis.

- ramr.method:

  A character scalar: when ramr.method is "IQR" (the default), the
  filtering based on interquantile range is used (\`iqr.cutoff\` value
  is then used as a threshold). When "beta", "wbeta" or "beinf" -
  filtering based on fitting non-weighted (\`EnvStats::ebeta\`),
  weighted (\`ExtDist::eBeta\`) or zero-and-one inflated
  (\`gamlss.dist::BEINF\`) beta distributions, respectively, is used,
  and \`pval.cutoff\` or \`qval.cutoff\` (if not \`NULL\`) is used as a
  threshold. For "wbeta", weights directly correlate with bin contents
  (number of values per bin) and inversly - with the distances from the
  median value, thus narrowing the estimated distribution and
  emphasizing outliers.

- iqr.cutoff:

  A single integer \>= 1. Methylation beta values differing from the
  median value by more than \`iqr.cutoff\` interquartile ranges are
  considered to be significant (the default: 5).

- pval.cutoff:

  A numeric scalar (the default: 5e-2). Bonferroni correction of
  \`pval.cutoff\` by the length of the \`data.samples\` object is used
  to calculate \`qval.cutoff\` if the latter is \`NULL\`.

- qval.cutoff:

  A numeric scalar. Used as a threshold for filtering based on fitting
  non-weighted or weighted beta distributions: all p-values lower than
  \`qval.cutoff\` are considered to be significant. If \`NULL\` (the
  default), it is calculated using \`pval.cutoff\`

- merge.window:

  A positive integer. All significant (survived the filtering stage)
  \`data.ranges\` genomic locations within this distance will be merged
  to create AMRs (the default: 300).

- min.cpgs:

  A single integer \>= 1. All AMRs containing less than \`min.cpgs\`
  significant genomic locations are filtered out (the default: 7).

- min.width:

  A single integer \>= 1 (the default). Only AMRs with the width of at
  least \`min.width\` are returned.

- exclude.range:

  A numeric vector of length two. If not \`NULL\` (the default), all
  \`data.ranges\` genomic locations with their median methylation beta
  value within the \`exclude.range\` interval are filtered out.

- cores:

  A single integer \>= 1. Number of processes for parallel computation
  (the default: all but one cores). Results of parallel processing are
  fully reproducible when the same seed is used (thanks to doRNG).

- verbose:

  boolean to report progress and timings (default: TRUE).

- ...:

  Further arguments to be passed to \`EnvStats::ebeta\` or
  \`ExtDist::eBeta\` functions.

## Value

The output is a \`GRanges\` object that contains all the aberrantly
methylated regions (AMRs) for all \`data.samples\` samples in
\`data.ranges\` object. The following metadata columns may be present:

- \`revmap\` – integer list of significant CpGs (\`data.ranges\` genomic
  locations) that are included in this AMR region

- \`ncpg\` – number of significant CpGs within this AMR region

- \`sample\` – contains an identifier of a sample to which corresponding
  AMR belongs

- \`dbeta\` – average deviation of beta values for significant CpGs from
  their corresponding median values

- \`pval\` – geometric mean of p-values for significant CpGs

- \`xiqr\` – average IQR-normalised deviation of beta values for
  significant CpGs from their corresponding median values

## Details

In the provided data set, \`getAMR.obsolete\` compares methylation beta
values of each sample with other samples to identify rare long-range
methylation aberrations (epimutations). For \`ramr.method=="IQR"\`: for
every genomic location (CpG) in \`data.ranges\` the IQR-normalized
deviation from the median value is calculated, and all CpGs with such
normalized deviation not smaller than the \`iqr.cutoff\` are retained.
For \`ramr.method distribution are estimated by means of
\`EnvStats::ebeta\` (beta distribution), \`ExtDist::eBeta\` (weighted
beta destribution), or \`gamlss.dist::BEINF\` (zero and one inflated
beta distribution) functions, respectively. These parameters are then
used to calculate the probability values, followed by the filtering when
all CpGs with p-values not greater than \`qval.cutoff\` are retained.
Another filtering is then performed to exclude all CpGs within
\`exclude.range\`. Next, the retained (significant) CpGs are merged
within the window of \`merge.window\`, and final filtering is applied to
AMR genomic ranges (by \`min.cpgs\` and \`min.width\`).

## See also

[`plotAMR`](plotAMR.md) for plotting AMRs,
[`getUniverse`](getUniverse.md) for info on enrichment analysis,
[`simulateAMR`](simulateAMR.md) and [`simulateData`](simulateData.md)
for the generation of simulated test data sets, and \`ramr\` vignettes
for the description of usage and sample data.

## Examples

``` r
  data(ramr)
  getAMR.obsolete(ramr.data, ramr.samples, ramr.method="beta",
                  min.cpgs=5, merge.window=1000, qval.cutoff=1e-3, cores=2)
#> Identifying AMRs
#> Loading required package: foreach
#> Loading required package: rngtools
#>  [5.423s]
#> Loading required namespace: GenomeInfoDb
#> GRanges object with 22 ranges and 5 metadata columns:
#>        seqnames          ranges strand |             revmap      ncpg
#>           <Rle>       <IRanges>  <Rle> |             <list> <integer>
#>    [1]     chr1 2443577-2453006      * | 2722,2723,2724,...        30
#>    [2]     chr1 1589891-1590941      * | 1459,1460,1461,...        10
#>    [3]     chr1 1589891-1590941      * | 1459,1460,1461,...        10
#>    [4]     chr1 1589891-1590941      * | 1459,1460,1461,...        10
#>    [5]     chr1   874697-877876      * |    165,166,167,...        13
#>    ...      ...             ...    ... .                ...       ...
#>   [18]     chr1 1709203-1715039      * | 1595,1596,1597,...        20
#>   [19]     chr1 1709203-1715039      * | 1595,1596,1597,...        20
#>   [20]     chr1   566172-569687      * |       17,18,19,...        15
#>   [21]     chr1 1138931-1146903      * |    726,727,728,...        27
#>   [22]     chr1 1160713-1165393      * |    789,790,791,...        15
#>             sample     dbeta         pval
#>        <character> <numeric>    <numeric>
#>    [1]    sample25 -0.484617  1.26772e-19
#>    [2]    sample33  0.489562  1.27388e-07
#>    [3]    sample34  0.503785  5.93348e-08
#>    [4]    sample35  0.512633  3.30858e-08
#>    [5]    sample44  0.451475  1.30458e-07
#>    ...         ...       ...          ...
#>   [18]    sample79  0.455621  2.80743e-23
#>   [19]    sample80  0.453449  3.24593e-23
#>   [20]    sample95  0.498337 4.60616e-171
#>   [21]    sample98  0.380622 3.95738e-232
#>   [22]    sample98 -0.514668  4.27297e-18
#>   -------
#>   seqinfo: 24 sequences from hg19 genome; no seqlengths
```
