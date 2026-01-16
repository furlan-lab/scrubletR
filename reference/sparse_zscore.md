# Z-score normalize each column of a sparse matrix

This function z-score normalizes each column of a sparse matrix.

## Usage

``` r
sparse_zscore(E, gene_mean = NULL, gene_stdev = NULL)
```

## Arguments

- E:

  A sparse matrix where rows represent observations and columns
  represent features.

- gene_mean:

  A vector of mean values for each feature (default is calculated from
  the matrix).

- gene_stdev:

  A vector of standard deviation values for each feature (default is
  calculated from the matrix).

## Value

A z-score normalized sparse matrix.
