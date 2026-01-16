# Calculate v-scores and related statistics for genes in the input counts matrix

This function calculates v-scores, coefficient of variation (CV), and
other statistics for genes in the input counts matrix. The v-score is an
above-Poisson noise statistic that helps assess the variability of gene
expression levels.

## Usage

``` r
get_vscores(E, min_mean = 0, nBins = 50, fit_percentile = 0.1, error_wt = 1)
```

## Arguments

- E:

  A counts matrix where rows represent cells and columns represent
  genes.

- min_mean:

  Minimum mean expression value for genes to be considered.

- nBins:

  Number of bins for calculating running quantiles.

- fit_percentile:

  Percentile used for fitting the running quantile.

- error_wt:

  Weight for the error function during optimization.

## Value

A list containing the following components:

- `v_scores`: V-scores for each gene.

- `CV_eff`: Coefficient of variation (CV) for effective counts.

- `CV_input`: Coefficient of variation (CV) for input counts.

- `gene_ix`: Indices of genes considered in the analysis.

- `mu_gene`: Mean expression values for selected genes.

- `FF_gene`: Fano factor (variance to mean ratio) for selected genes.

- `a`: Parameter 'a' obtained during optimization.

- `b`: Parameter 'b' obtained during optimization.

## See also

[`runningquantile`](https://furlan-lab.github.io/scrubletR/reference/runningquantile.md),
[`optimize`](https://rdrr.io/r/stats/optimize.html)
