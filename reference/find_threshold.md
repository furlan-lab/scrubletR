# Find the threshold between two maxima in a smoothed histogram

This function finds the threshold between two maxima in a smoothed
histogram. It iteratively smooths the histogram and then identifies the
maxima, breaking if fewer than three maxima are found.

## Usage

``` r
find_threshold(input, nbins = 256, max_num_iter = 100)
```

## Arguments

- nbins:

  Number of bins for histogram

- max_num_iter:

  Maximum number of iterations for histogram smoothing

- image:

  Input image

- hist:

  Numeric vector representing the histogram

## Value

The threshold between two maxima in the smoothed histogram
