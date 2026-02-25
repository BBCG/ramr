# Simulated Illumina HumanMethylation 450k data set with 3000 CpGs and 100 samples

Data was simulated using GSE51032 data set as described in the
reference. Current data set (`"ramr.data"`) contains beta values for
10000 CpGs and 100 samples (`"ramr.samples"`), and carries 6 unique
(`"ramr.tp.unique"`) and 15 non-unique (`"ramr.tp.nonunique"`) true
positive AMRs containing at least 10 CpGs with their beta values
increased/decreased by 0.5.

## Usage

``` r
data(ramr)
```

## Format

Objects of class `"GRanges"`
(`"ramr.data, ramr.tp.unique, ramr.tp.nonunique"`) and `"character"`
(`"ramr.samples"`).

## References

Nikolaienko et al., 2020
([bioRxiv](https://doi.org/10.1101/2020.12.01.403501))

## Examples

``` r
  data(ramr)
  amrs <- getAMR(
    data.ranges=ramr.data, compute="IQR",
    combine.min.cpgs=5, combine.window=1000, combine.threshold=5
  )
#> Preprocessing data 
#> [0.056s]
#> Computing IQR 
#> [0.003s]
#> Creating genomic ranges 
#> [0.009s]
  plotAMR(data.ranges=ramr.data, amr.ranges=amrs[1])
#> Plotting 1 genomic ranges 
#> 100%
#> [0.160s]
#> $`chr1:2443577-2453006`

#> 
  plotAMR(data.ranges=ramr.data, amr.ranges=ramr.tp.nonunique[4],
          highlight=c("sample7","sample8","sample9"))
#> Plotting 1 genomic ranges 
#> 100%
#> [0.159s]
#> $`chr1:2119531-2122601`

#> 
```
