# Subsample counts in a sparse matrix

This function subsamples counts in a sparse matrix based on a given
rate.

## Usage

``` r
subsample_counts(E, rate, original_totals, random_seed = 0)
```

## Arguments

- E:

  A sparse matrix where rows represent observations and columns
  represent features.

- rate:

  Subsampling rate for counts (values between 0 and 1).

- original_totals:

  A vector of original total counts for each observation.

- random_seed:

  An integer specifying the random seed for reproducibility.

## Value

A list containing the subsampled sparse matrix (E) and the final
downsampling totals.
