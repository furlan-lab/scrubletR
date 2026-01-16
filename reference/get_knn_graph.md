# Get k-Nearest-Neighbor Graph

Build k-nearest-neighbor graph and return edge list and nearest neighbor
matrix.

## Usage

``` r
get_knn_graph(
  X,
  k = 5,
  dist_metric = c("euclidean", "angular", "manhattan", "hamming"),
  approx = FALSE,
  return_edges = TRUE,
  random_seed = 0
)
```

## Arguments

- X:

  Data matrix.

- k:

  Number of nearest neighbors.

- dist_metric:

  Distance metric for finding neighbors.

- approx:

  Whether to use approximate nearest neighbor search.

- return_edges:

  Whether to return edge list.

- random_seed:

  Random seed for reproducibility.

## Value

If `return_edges` is `TRUE`, a list containing edge list and nearest
neighbor matrix; otherwise, the nearest neighbor matrix.
