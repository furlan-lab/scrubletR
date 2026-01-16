# Scrublet

Wolock SL, Lopez R, Klein AM. Scrublet: Computational Identification of
Cell Doublets in Single-Cell Transcriptomic Data. Cell Syst. 2019 Apr
24;8(4):281-291.e9. doi: 10.1016/j.cels.2018.11.005. Epub 2019 Apr 3.
PMID: 30954476; PMCID: PMC6625319.
https://www.sciencedirect.com/science/article/pii/S2405471218304745

## Usage

``` r
scrublet(
  object,
  split_by = NULL,
  return_results_only = FALSE,
  min_counts = 3,
  min_cells = 3,
  min_gene_variability_pctl = 85,
  seed = 2024,
  expected_doublet_rate = 0.1,
  n_prin_comps = 30,
  sim_doublet_ratio = 2,
  assay = "RNA",
  cores = 1,
  show_gene_filter_plot = F
)
```

## Arguments

- object:

  the object upon which to perform Scrublet (monocle3 objects and seurat
  supported)

- split_by:

  the column in the meta data to split the object by before running
  scrublet

- return_results_only:

  bool (optional, default False)

- min_counts, :

  int (optional, default=2), See scrublet reference

- min_cells, :

  int (optional, default=3), See scrublet reference

- min_gene_variability_pctl, :

  int (optional, default=85), See scrublet reference

- seed, :

  seed aka random state

- expected_doublet_rate:

  = float (default 0.1); doesn't affect doublet score calculation only
  prediction.

- n_prin_comps, :

  int (optional, default=30), See scrublet reference (Number of
  principal components to use)

- sim_doublet_ratio, :

  int (optional, default=2), the number of doublets to simulate,
  relative to the number of observed transcriptomes. This should be high
  enough that all doublet states are well-represented by simulated
  doublets. Setting it too high is computationally expensive. The
  default value is 2, though values as low as 0.5 give very similar
  results for the datasets that have been tested.

- cores:

  Number of cores (only helps when splitting and object)

- show_gene_filter_plot, :

  show the mean x FF plot with selected features

## Value

The input CellDataSet or Seurat object with an additional column added
to pData with both the doublet_score output from scrublet,
