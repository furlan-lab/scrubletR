# Total Counts Normalization

Total Counts Normalization

## Usage

``` r
tot_counts_norm(
  E,
  total_counts = NULL,
  exclude_dominant_frac = 1,
  included = integer(0),
  target_total = NULL
)
```

## Arguments

- E:

  Counts matrix to be normalized

- total_counts:

  Vector of total counts per cell (if NULL, computed from E)

- exclude_dominant_frac:

  Exclude overly abundant genes fraction (default: 1)

- included:

  Indices of genes to include in normalization (if empty, all genes are
  considered)

- target_total:

  Target total for normalization (if NULL, the mean of total counts is
  used)

## Value

List containing normalized counts matrix (Enorm), average total counts,
and vector of included genes
