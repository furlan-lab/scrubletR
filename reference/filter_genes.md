# Filter genes by expression level and variability

This function filters genes based on expression level and variability
using v-scores.

## Usage

``` r
filter_genes(
  E,
  base_ix = NULL,
  min_vscore_pctl = 85,
  min_counts = 3,
  min_cells = 3,
  plot = FALSE,
  sample_name = ""
)
```

## Arguments

- E:

  A counts matrix where rows represent cells and columns represent
  genes.

- base_ix:

  Indices of cells to be used for v-score calculation (default: all
  cells).

- min_vscore_pctl:

  Minimum percentile threshold for v-scores.

- min_counts:

  Minimum expression counts required for a gene to be considered.

- min_cells:

  Minimum number of cells expressing a gene for it to be considered.

- sample_name:

  Character string specifying the name of the sample for the plot title.

- show_vscore_plot:

  Logical indicating whether to show a v-score plot.

## Value

A numeric vector containing the indices of filtered genes.

## See also

[`get_vscores`](https://furlan-lab.github.io/scrubletR/reference/get_vscores.md)
