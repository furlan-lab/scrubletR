# Find local maxima indices in a histogram

This function finds local maxima indices in a histogram, as scipy's
argrelmax fails on plateaus.

## Usage

``` r
find_local_maxima_idx(hist)
```

## Arguments

- hist:

  Numeric vector representing the histogram

## Value

A list containing indices of local maxima in the histogram
