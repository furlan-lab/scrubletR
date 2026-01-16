# Calculate running quantiles for given data points

This function calculates running quantiles for given x and y data
points.

## Usage

``` r
runningquantile(x, y, p, nBins)
```

## Arguments

- x:

  Numeric vector representing the x-axis data points.

- y:

  Numeric vector representing the y-axis data points.

- p:

  Percentile value for calculating the running quantiles.

- nBins:

  Number of bins for the running quantiles.

## Value

A list containing the following components:

- `xOut`: x-axis values for the running quantiles.

- `yOut`: y-axis values representing the running quantiles.
