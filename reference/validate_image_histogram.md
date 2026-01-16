# Helper function to validate image histogram.

Helper function to validate image histogram.

## Usage

``` r
validate_image_histogram(image, hist, nbins)
```

## Arguments

- image:

  numeric vector Grayscale input image.

- hist:

  array, or 2-tuple of arrays Histogram to determine the threshold from
  and a corresponding array of bin center intensities.

- nbins:

  int Number of bins used to calculate histogram.

## Value

list A list containing counts and bin_centers.
