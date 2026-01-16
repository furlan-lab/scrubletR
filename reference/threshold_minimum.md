# Return threshold value based on minimum method.

The histogram of the input `image` is computed if not provided and
smoothed until there are only two maxima. Then the minimum in between is
the threshold value.

## Usage

``` r
threshold_minimum(image = NULL, nbins = 256, max_num_iter = 10000, hist = NULL)
```

## Arguments

- image:

  (M, N\[, ...\]) numeric vector, optional Grayscale input image.

- nbins:

  int, optional Number of bins used to calculate histogram. This value
  is ignored for integer arrays.

- max_num_iter:

  int, optional Maximum number of iterations to smooth the histogram.

- hist:

  array, or 2-tuple of arrays, optional Histogram to determine the
  threshold from and a corresponding array of bin center intensities.
  Alternatively, only the histogram can be passed.

## Value

threshold float Upper threshold value. All pixels with an intensity
higher than this value are assumed to be foreground.

## Details

Either image or hist must be provided. In case hist is given, the actual
histogram of the image is ignored.

## References

C. A. Glasbey, "An analysis of histogram-based thresholding algorithms,"
CVGIP: Graphical Models and Image Processing, vol. 55, pp. 532-537,
1993. Prewitt, JMS & Mendelsohn, ML (1966), "The analysis of cell
images", Annals of the New York Academy of Sciences 128: 1035-1053
:DOI:\`10.1111/j.1749-6632.1965.tb11715.x\`
