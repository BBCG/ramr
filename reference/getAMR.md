# Search for aberrantly methylated regions

\`getAMR\` returns a \`GRanges\` object with aberrantly methylated
regions (AMRs / epimutations) for all samples in a data set.

## Usage

``` r
getAMR(
  data.ranges,
  data.samples = NULL,
  data.coverage = NULL,
  transform = c("identity", "linear"),
  exclude.range = NULL,
  compute = c("IQR", "beta+binom"),
  compute.estimate = c("mom", "amle", "nmle"),
  compute.weights = c("equal", "logInvDist", "sqrtInvDist", "invDist"),
  combine = c("threshold", "comb-p"),
  combine.threshold = ifelse(compute == "IQR", 5, 0.001),
  combine.window = 300,
  combine.min.cpgs = 7,
  combine.min.width = 1,
  combine.ignore.strand = FALSE,
  ncores = NULL,
  verbose = TRUE
)
```

## Arguments

- data.ranges:

  A \`GRanges\` object with genomic locations and corresponding beta
  values included as metadata.

- data.samples:

  A character vector with sample names (e.g., a subset of metadata
  column names). If \`NULL\` (the default), then all samples (metadata
  columns) are included in the analysis.

- data.coverage:

  An optional \`data.frame\` object with coverage data. If provided,
  must be an all-integer \`data.frame\`, with the same dimensions and
  column names as the \`data.ranges\` metadata columns.

- transform:

  A character scalar specifying if beta values should be used as
  supplied ("identity", the default) or linearly transformed ("linear")
  to force all \\\\0;1\\\\ endpoint (extreme) values within open
  \\(0,1)\\ interval using the following formula:
  \$\$x'=\frac{x(N-1)+0.5}{N}\$\$ where \\x\\ is the beta value before
  transformation, and \\N\\ is the number of samples. Such
  transformation is recommended only if beta values contain \\\\0;1\\\\
  values **and** coverage data is NOT available. See [doi:
  10.1037/1082-989X.11.1.54](https://doi.org/10.1037/1082-989x.11.1.54)
  and [Dealing with 0,1 values in a beta
  regression](https://stats.stackexchange.com/questions/31300/dealing-with-0-1-values-in-a-beta-regression)
  for more details.

- exclude.range:

  A numeric vector of length two. Unless \`NULL\` (the default), all
  \`data.ranges\` genomic locations with their median methylation beta
  value within the \`exclude.range\` interval are filtered out.

- compute:

  A character scalar: when "IQR" (the default), AMR search based on
  interquantile range is performed. When "beta+binom" - search is based
  on fitting optionally weighted beta distribution (with the possibility
  of estimating binomial probability of \\\\0;1\\\\ endpoint values).
  See Details section for explanation.

- compute.estimate:

  A character scalar for the method of parameter estimation of beta
  distribution. The default ("mom") stands for the method of moments
  based on the unbiased estimator of variance and includes \\\\0;1\\\\
  endpoints in calculation of moments (mean, unbiased variance). Other
  options are "amle" (approximation of maximum likelihood estimation)
  and "nmle" (numeric maximum likelihood estimation) - both ignore
  \\\\0;1\\\\ endpoints in calculations. See Details section for
  explanation.

- compute.weights:

  A character scalar for the weights assigned to individual observations
  (beta values) and used to compute weighted means and variance as
  described in "Compute" section below. Four available weighing schemes
  result in different sensitivity of outlier detection and rate of false
  positive (FP) findings: "equal" (the default) is the least sensitive
  and gives least number of FPs, while "invDist" is the most sensitive
  but may result in a very high number of FPs especially when
  \`combine.threshold\` is too high (1e-3 or higher). "logInvDist" is
  recommended when one desires a balance between relatively low type I
  error rate and higher detection sensitivity for both unique and
  non-unique AMRs.

  If "equal", all weights are equal to 1 (\\w_i=1\\). Otherwise, weights
  of observations inversely depend on their distance from the median
  (thus emphasizing outliers) and are calculated using the following
  formulas: \$\$\text{"logInvDist": } w_i =
  \log{\frac{1}{\|M-x_i\|+\epsilon}}\$\$ \$\$\text{"sqrtInvDist": } w_i
  = \sqrt{\frac{1}{\|M-x_i\|+\epsilon}}\$\$ \$\$\text{"invDist": } w_i =
  \frac{1}{\|M-x_i\|+\epsilon}\$\$ where \\x_i\\ is a beta value for
  i-th sample, \\M\\ is a median of all beta values at this genomic
  location, and \\\epsilon\\ is a very small number (= FLT_EPSILON
  \\\approx\\ 1.192093e-07).

- combine:

  A character scalar for the method used to combine individual outlier
  genomic positions into genomic intervals (ranges). When "threshold"
  (the default) is used, simple thresholding is applied. More details
  are given in the "Combine" subsection below.

- combine.threshold:

  A numeric scalar setting the threshold for an outlier value. When
  \`compute=="IQR"\`, methylation beta values differing from the median
  value by at least \`combine.threshold\` interquartile ranges are
  considered to be outliers (the default: 5). When
  \`compute=="beta+binom"\`, all probability values not higher than
  \`combine.threshold\` are considered to be outliers (the default:
  0.001).

- combine.window:

  A positive integer. All significant (survived the filtering stage)
  \`data.ranges\` genomic locations within this distance will be merged
  to create AMRs (the default: 300).

- combine.min.cpgs:

  A single integer \>= 1. All AMRs containing less than
  \`combine.min.cpgs\` significant genomic locations are filtered out
  (the default: 7).

- combine.min.width:

  A single integer \>= 1 (the default). Only AMRs with the width of at
  least \`combine.min.width\` are returned.

- combine.ignore.strand:

  A boolean scalar to ignore strand information in the input
  \`data.ranges\` and when outlier genomic positions are combined into
  genomic intervals (the default: FALSE).

- ncores:

  A single integer \>= 1. Number of OpenMP threads for parallel
  computation (the default: half of the available cores). Results of
  parallel processing are fully reproducible.

- verbose:

  Boolean to report progress and timings (default: TRUE).

## Value

The output is a \`GRanges\` object that contains all the aberrantly
methylated regions (AMRs / epimutations) for all \`data.samples\`
samples in \`data.ranges\` object. The following metadata columns may be
present:

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

In the provided data set, \`getAMR\` finds stretches of outlier beta
values to identify rare long-range methylation aberrations
(epimutations) in one or several samples. As a rule, methods for
differential methylation analysis rely on between-group comparisons —
\`getAMR\` performs this comparison within-group, which is not only
faster, but also more sensitive. The logic of computations is described
below.

### Compute

This section describes computations that are performed in order to
identify individual outlier beta values — but **before values are deemed
as outliers**. At the moment, the two supported methods are:

- "IQR":

  When \`compute=="IQR"\`, for every genomic location (CpG) in
  \`data.ranges\` the IQR-normalized deviation from the median value is
  calculated using the following formula:
  \$\$xIQR_i=\frac{x_i-M}{IQR}\$\$ where \\x_i\\ is a beta value for
  i-th sample, \\M\\ and \\IQR\\ are a median and an interquartile range
  of all beta values at this genomic location, respectively.

- "beta+binom":

  When \`compute=="beta+binom"\`, for every genomic location (CpG),
  \`getAMR\` will estimate the probability of each beta value to occur.
  For all beta values inside the open \\(0,1)\\ interval, beta
  distribution is used. For all \\\\0;1\\\\ endpoint (extreme) values
  where beta distribution is not defined, binomial probability is
  calculated using coverage data (supplied using \`data.coverage\`
  parameter).

  Alpha \\{\alpha}\\ and beta \\{\beta}\\ parameters of beta
  distribution can be estimated using one of the following methods:

  1.  Method of moments (\`compute.estimate="mom"\`) based on (both
      optionally weighted) mean and unbiased variance \$\$\text{sample
      mean} = \bar{x} = \frac{ \sum\limits\_{i=1}^n w_i
      x_i}{\sum\limits\_{i=1}^n w_i}\$\$ \$\$\text{unbiased variance}
      =\bar{v} = \frac {\sum\limits\_{i=1}^N w_i (x_i - \bar{x})^2}
      {\sum\_{i=1}^N w_i - (\sum\_{i=1}^N w_i^2 / \sum\_{i=1}^N w_i)}
      \$\$ where \\x_i\\ and \\w_i\\ are a beta value and its
      reliability weight for i-th sample. \$\$\hat{\alpha} = \bar{x}
      \left(\frac{\bar{x} (1 - \bar{x})}{\bar{v}} - 1 \right)\$\$
      \$\$\hat{\beta} = (1-\bar{x}) \left(\frac{\bar{x} (1 -
      \bar{x})}{\bar{v}} - 1 \right)\$\$ Note: all beta values from the
      closed \\\[0,1\]\\ interval are used for calculation of arithmetic
      mean value. For more details on used formulas, visit [Weighted
      arithmetic
      mean](https://en.wikipedia.org/wiki/Weighted_arithmetic_mean#Mathematical_definition),
      [Weighted variance with reliability
      weights](https://en.wikipedia.org/wiki/Weighted_arithmetic_mean#Reliability_weights),
      and [Method of moments for beta
      distribution](https://en.wikipedia.org/wiki/Beta_distribution#Method_of_moments)

  2.  Approximate maximum likelihood (\`compute.estimate="amle"\`) based
      on (both optionally weighted) geometric means of \\x\\ and
      \\(1-x)\\ \$\$ \hat{G}\_x = \left(\prod\_{i=1}^n
      x_i^{w_i}\right)^{1 / \sum\_{i=1}^n w_i} \$\$ \$\$
      \hat{G}\_{(1-x)} = \left(\prod\_{i=1}^n {(1-x_i)}^{w_i}\right)^{1
      / \sum\_{i=1}^n w_i} \$\$ where \\x_i\\ and \\w_i\\ are a beta
      value and its reliability weight for i-th sample. \$\$
      \hat{\alpha}\approx \tfrac{1}{2} +
      \frac{\hat{G}\_{x}}{2(1-\hat{G}\_x-\hat{G}\_{(1-x)})} \$\$ \$\$
      \hat{\beta}\approx \tfrac{1}{2} +
      \frac{\hat{G}\_{(1-x)}}{2(1-\hat{G}\_x-\hat{G}\_{(1-x)})} \$\$
      Note: only beta values from the open \\(0,1)\\ interval are used
      for calculation of geometric means. For more details, visit
      [Weighted geometric
      mean](https://en.wikipedia.org/wiki/Weighted_geometric_mean), and
      [Maximum likelihood for beta
      distribution](https://en.wikipedia.org/wiki/Beta_distribution#Two_unknown_parameters_2)

  3.  And a numerical maximum likelihood (\`compute.estimate="nmle"\`)
      which is yet to be implemented.

  After estimating parameters of beta distribution, probabilities of
  observed beta values from the open \\(0,1)\\ interval are computed
  using regularized incomplete beta function \\I_x(\alpha, \beta)\\. For
  more details, visit [Cumulative distribution function for beta
  distribution](https://en.wikipedia.org/wiki/Beta_distribution#Cumulative_distribution_function).

  Probabilities of \\\\0;1\\\\ endpoint values are computed using mean
  value (either arithmetic for \`compute.estimate="mom"\` or geometric
  for \`compute.estimate="\*mle"\`) using the following formulas:
  \$\${p^0}\_i=(1-\bar{x})^k\$\$ \$\${p^1}\_i=\bar{x}^k\$\$ where
  \\{p^0}\_i\\ and \\{p^1}\_i\\ are probabilities of observing 0 or 1
  for i-th sample, respectively, \\\bar{x}\\ is a mean of all beta
  values for this genomic position, and \\k\\ is a sequencing coverage
  of this genomic position for i-th sample.

### Combine

This section describes how individual beta values are deemed as outliers
and how these outliers are combined into extended genomic regions
(aberrantly methylated regions, AMRs, or epimutations). The only method
supported at the moment is thresholding as described below.

- "threshold":

  This method applies a fixed threshold (specified using
  \`combine.threshold\` parameter) to all \`xIQR\` or probability values
  computed at the previous step. All observed beta values that are
  passing this threshold (above or equal to the threshold for \`xIQR\`;
  below or equal in case of probability values) are considered to be
  outliers and therefore retained.

  Next, for all outlier beta values per sample, corresponding genomic
  positions are merged into genomic intervals using the window of
  \`combine.window\` and keeping the strand information unless
  \`combine.ignore.strand\` is TRUE.

- "comb-p":

  This method of combining outliers into genomic ranges is yet to be
  implemented.

Resulting genomic intervals are then filtered: only the regions
containing at least \`combine.min.cpgs\` outliers and which are at leas
as wide as \`combine.min.width\` are reported back.

As a final step, the following average values are computed for every
aberrantly methylated region:

1.  arithmetic mean of distances of all outlier beta values to a median
    beta value (\`dbeta\`)

2.  if \`compute=="IQR"\`, arithmetic mean of xIQR values of all outlier
    genomic position (\`xiqr\`) OR

3.  if \`compute=="beta+binom"\`, geometric mean of probability values
    of all outlier genomic position (\`pval\`)

## Note

NA values within metadata columns of \`data.ranges\` are silently
dropped in all computations.

## See also

[`plotAMR`](plotAMR.md) for plotting AMRs,
[`getUniverse`](getUniverse.md) for info on enrichment analysis,
[`simulateAMR`](simulateAMR.md) and [`simulateData`](simulateData.md)
for the generation of simulated test data sets, and \`ramr\` vignettes
for the description of usage and sample data.

## Examples

``` r
  data(ramr)
  getAMR(data.ranges=ramr.data, data.samples=ramr.samples,
         compute="beta+binom", compute.estimate="amle",
         combine.min.cpgs=5, combine.window=1000, combine.threshold=1e-3)
#> Preprocessing data 
#> [0.048s]
#> Fitting beta distribution 
#> [0.041s]
#> Creating genomic ranges 
#> [0.022s]
#> GRanges object with 21 ranges and 5 metadata columns:
#>        seqnames          ranges strand |             revmap      ncpg   sample
#>           <Rle>       <IRanges>  <Rle> |             <list> <integer> <factor>
#>    [1]     chr1 2443577-2453006      * | 2722,2723,2724,...        30 sample25
#>    [2]     chr1 1589891-1590941      * | 1459,1460,1461,...        10 sample33
#>    [3]     chr1 1589891-1590941      * | 1459,1460,1461,...        10 sample34
#>    [4]     chr1 1589891-1590941      * | 1459,1460,1461,...        10 sample35
#>    [5]     chr1   874697-877876      * |    165,166,167,...        13 sample44
#>    ...      ...             ...    ... .                ...       ...      ...
#>   [17]     chr1 1709203-1715039      * | 1595,1596,1597,...        20 sample79
#>   [18]     chr1 1709203-1715039      * | 1595,1596,1597,...        20 sample80
#>   [19]     chr1   566172-569687      * |       17,18,19,...        15 sample95
#>   [20]     chr1 1138931-1146903      * |    726,727,728,...        27 sample98
#>   [21]     chr1 1160713-1165393      * |    789,790,791,...        15 sample98
#>            dbeta        pval
#>        <numeric>   <numeric>
#>    [1] -0.484617 1.13977e-19
#>    [2]  0.489562 8.39634e-08
#>    [3]  0.503785 3.82302e-08
#>    [4]  0.512633 2.10071e-08
#>    [5]  0.451475 9.63767e-08
#>    ...       ...         ...
#>   [17]  0.455621 1.76597e-08
#>   [18]  0.453449 2.04956e-08
#>   [19]  0.498337 6.26893e-17
#>   [20]  0.380622 1.77707e-18
#>   [21] -0.514668 2.65810e-18
#>   -------
#>   seqinfo: 24 sequences from an unspecified genome; no seqlengths
```
