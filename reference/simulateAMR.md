# Simulate a set of aberrantly methylated regions

\`simulateAMR\` returns a \`GRanges\` object containing a set of
randomly selected aberrantly methylated regions (AMRs) to be used as an
input for the \`simulateData\` method.

## Usage

``` r
simulateAMR(
  template.ranges,
  nsamples,
  exclude.ranges = NULL,
  regions.per.sample = 1,
  samples.per.region = 1,
  sample.names = NULL,
  merge.window = 300,
  min.cpgs = 7,
  max.cpgs = Inf,
  min.width = 1,
  dbeta = 0.25
)
```

## Arguments

- template.ranges:

  A \`GRanges\` object with genomic locations (same object must be
  supplied to this and to the \`simulateData\` functions).

- nsamples:

  A single integer \>= 1 indicating the number of samples to which AMRs
  will be assigned.

- exclude.ranges:

  A \`GRanges\` object with genomic locations. None of the simulated
  AMRs in the output will overlap with any of regions from
  \`exclude.ranges\`. If \`NULL\` (the default), AMRs are not restricted
  by their genomic location.

- regions.per.sample:

  A single integer \>= 1 (the default). Number of AMRs to be assigned to
  every sample. Message is shown and the \`regions.per.sample\` value is
  limited to \`maxnAMR (where \`maxnAMR\` is the maximum number of
  potential AMRs for the \`template.ranges\`).

- samples.per.region:

  A single integer \>= 1 (the default). Number of samples to which the
  same AMR will be assigned. Message is shown and the
  \`samples.per.region\` value is limited to \`nsamples\` if the former
  is greater than the latter.

- sample.names:

  A character vector with sample names. If \`NULL\` (the default),
  sample names will be computed as \`paste0("sample",
  seq_len(nsamples))\`. When specified, the length of the
  \`sample.names\` vector must not be smaller than the value of
  \`nsamples\`.

- merge.window:

  A positive integer. All \`template.ranges\` genomic locations within
  this distance will be merged to create a list of potential AMRs (which
  will be later filtered from regions overlapping with any regions from
  the \`exclude.ranges\`).

- min.cpgs:

  A single integer \>= 1. All AMRs containing less than \`min.cpgs\`
  genomic locations are filtered out. The default: 7.

- max.cpgs:

  A single integer \>= 1. All AMRs containing more than \`max.cpgs\`
  genomic locations are filtered out. The default: \`Inf\`.

- min.width:

  A single integer \>= 1 (the default). Only AMRs with the width of at
  least \`min.width\` are returned.

- dbeta:

  A single non-negative numeric value in the range \[0,1\] or a numeric
  vector of such values (with as many elements as there are AMRs). Used
  to populate the \`dbeta\` metadata column, defines a desired absolute
  deviation of corresponding AMR from the median for the
  \`simulateData\` function.

## Value

The output is a \`GRanges\` object that contains a subset of aberrantly
methylated regions (AMRs) randomly selected from all the possible AMRs
for the provided \`template.ranges\` object. The following metadata
columns are included:

- \`revmap\` – integer list of \`template.ranges\` genomic locations
  that are included in this AMR region

- \`ncpg\` – number of \`template.ranges\` genomic locations within this
  AMR region

- \`sample\` – an identifier of a sample to which corresponding AMR
  belongs

- \`dbeta\` – equals to supplied \`dbeta\` parameter

## Details

Using provided template (\`GRanges\` object) \`simulateAMR\` randomly
selects genomic regions satisfying various criteria (number of CpGs,
width of the region) and assigns them to samples according to specified
parameters (number of AMRs per sample, number of samples per AMR). Its
output is meant to be used as the set of true positive AMRs for the
\`simulateData\` function.

## See also

[`simulateData`](simulateData.md) for the generation of simulated test
data sets, [`getAMR`](getAMR.md) for identification of AMRs,
[`plotAMR`](plotAMR.md) for plotting AMRs,
[`getUniverse`](getUniverse.md) for info on enrichment analysis, and
\`ramr\` vignettes for the description of usage and sample data.

## Examples

``` r
  data(ramr)
  set.seed(1)
  amrs.unique <-
    simulateAMR(ramr.data, nsamples=4, regions.per.sample=2,
                min.cpgs=5, merge.window=1000, dbeta=0.2)
  amrs.nonunique <-
    simulateAMR(ramr.data, nsamples=3, exclude.ranges=amrs.unique,
                samples.per.region=2, min.cpgs=5, merge.window=1000)
```
