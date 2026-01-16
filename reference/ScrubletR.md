# Scrublet R6 Class

The Scrublet R6 class represents a single-cell RNA-seq scrubbing tool.

## Public Fields

- `E_obs`: Observed expression matrix.

- `E_sim`: Simulated expression matrix.

- `E_obs_norm`: Normalized observed expression matrix.

- `E_sim_norm`: Normalized simulated expression matrix.

- `gene_filter`: Filter for highly variable genes.

- `embeddings`: List of embeddings.

- `total_counts_obs`: Total counts for observed cells.

- `total_counts_sim`: Total counts for simulated cells.

- `sim_doublet_ratio`: Ratio of simulated doublets to observed cells.

- `n_neighbors`: Number of neighbors for calculations.

- `expected_doublet_rate`: Expected doublet rate.

- `stdev_doublet_rate`: Standard deviation of doublet rate.

- `random_state`: Seed for reproducibility.

- `doublet_parents`: Matrix of doublet parents.

- `manifold_obs_`: Observed manifold data.

- `manifold_sim_`: Simulated manifold data.

- `doublet_scores_obs_`: Doublet scores for observed cells.

- `doublet_scores_sim_`: Doublet scores for simulated cells.

- `doublet_errors_obs_`: Doublet errors for observed cells.

- `doublet_errors_sim_`: Doublet errors for simulated cells.

- `doublet_neighbor_parents_`: Parents of doublet neighbors.

- `predicted_doublets`: Predicted doublets.

- `z_scores_`: Z-scores for doublet predictions.

- `threshold_`: Threshold for doublet predictions.

- `detected_doublet_rate_`: Detected doublet rate.

- `detectable_doublet_fraction_`: Detectable doublet fraction.

- `overall_doublet_rate_`: Overall doublet rate.

## Public Methods

- `get_dims`: Display dimensions of expression matrices.

- `scrub_doublets`: Perform scrubbing to identify and remove doublets.

- `simulate_doublets`: Simulate doublets based on observed data.

- `set_manifold`: Set manifold data.

- `calculate_doublet_scores`: Calculate doublet scores using nearest
  neighbors.

- `call_doublets`: Identify doublets based on scores and threshold.

- `nearest_neighbor_classifier`: Nearest neighbor classification for
  doublet scores.

- `plot_histogram`: Plot histogram of doublet scores.

- `set_embedding`: Set embedding data.

- `plot_embedding`: Plot embeddings.

- `pipeline_normalize`: Normalize total counts.

- `pipeline_get_gene_filter`: Identify highly variable genes.

- `pipeline_apply_gene_filter`: Apply gene filter to expression
  matrices.

- `pipeline_mean_center`: Mean center expression matrix.

- `pipeline_normalize_variance`: Variance normalization of expression
  matrices.

- `pipeline_zscore`: Z-score normalization of expression matrices.

- `pipeline_log_transform`: Log transform expression matrices.

- `pipeline_truncated_svd`: Truncated Singular Value Decomposition.

- `pipeline_pca`: Principal Component Analysis.

## Methods

### Public methods

- [`ScrubletR$get_dims()`](#method-Scrublet-get_dims)

- [`ScrubletR$new()`](#method-Scrublet-new)

- [`ScrubletR$scrub_doublets()`](#method-Scrublet-scrub_doublets)

- [`ScrubletR$simulate_doublets()`](#method-Scrublet-simulate_doublets)

- [`ScrubletR$set_manifold()`](#method-Scrublet-set_manifold)

- [`ScrubletR$calculate_doublet_scores()`](#method-Scrublet-calculate_doublet_scores)

- [`ScrubletR$call_doublets()`](#method-Scrublet-call_doublets)

- [`ScrubletR$nearest_neighbor_classifier()`](#method-Scrublet-nearest_neighbor_classifier)

- [`ScrubletR$plot_histogram()`](#method-Scrublet-plot_histogram)

- [`ScrubletR$set_embedding()`](#method-Scrublet-set_embedding)

- [`ScrubletR$plot_embedding()`](#method-Scrublet-plot_embedding)

- [`ScrubletR$pipeline_normalize()`](#method-Scrublet-pipeline_normalize)

- [`ScrubletR$pipeline_get_gene_filter()`](#method-Scrublet-pipeline_get_gene_filter)

- [`ScrubletR$pipeline_apply_gene_filter()`](#method-Scrublet-pipeline_apply_gene_filter)

- [`ScrubletR$pipeline_mean_center()`](#method-Scrublet-pipeline_mean_center)

- [`ScrubletR$pipeline_normalize_variance()`](#method-Scrublet-pipeline_normalize_variance)

- [`ScrubletR$pipeline_zscore()`](#method-Scrublet-pipeline_zscore)

- [`ScrubletR$pipeline_log_transform()`](#method-Scrublet-pipeline_log_transform)

- [`ScrubletR$pipeline_truncated_svd()`](#method-Scrublet-pipeline_truncated_svd)

- [`ScrubletR$pipeline_pca()`](#method-Scrublet-pipeline_pca)

- [`ScrubletR$clone()`](#method-Scrublet-clone)

------------------------------------------------------------------------

### Method `get_dims()`

#### Usage

    ScrubletR$get_dims()

------------------------------------------------------------------------

### Method [`new()`](https://rdrr.io/r/methods/new.html)

#### Usage

    ScrubletR$new(
      counts_matrix,
      total_counts = NULL,
      sim_doublet_ratio = 2,
      n_neighbors = NULL,
      expected_doublet_rate = 0.1,
      stdev_doublet_rate = 0.02,
      random_state = 0,
      show_gene_filter_plot = TRUE
    )

------------------------------------------------------------------------

### Method `scrub_doublets()`

#### Usage

    ScrubletR$scrub_doublets(
      synthetic_doublet_umi_subsampling = 1,
      use_approx_neighbors = TRUE,
      distance_metric = "euclidean",
      get_doublet_neighbor_parents = FALSE,
      min_counts = 3,
      min_cells = 3,
      min_gene_variability_pctl = 85,
      log_transform = FALSE,
      mean_center = T,
      normalize_variance = T,
      n_prin_comps = 30,
      verbose = TRUE
    )

------------------------------------------------------------------------

### Method `simulate_doublets()`

#### Usage

    ScrubletR$simulate_doublets(
      sim_doublet_ratio = NULL,
      synthetic_doublet_umi_subsampling = 1
    )

------------------------------------------------------------------------

### Method `set_manifold()`

#### Usage

    ScrubletR$set_manifold(manifold_obs, manifold_sim)

------------------------------------------------------------------------

### Method `calculate_doublet_scores()`

#### Usage

    ScrubletR$calculate_doublet_scores(
      use_approx_neighbors = TRUE,
      distance_metric = "euclidean",
      get_doublet_neighbor_parents = FALSE
    )

------------------------------------------------------------------------

### Method `call_doublets()`

#### Usage

    ScrubletR$call_doublets(threshold = NULL, verbose = TRUE)

------------------------------------------------------------------------

### Method `nearest_neighbor_classifier()`

#### Usage

    ScrubletR$nearest_neighbor_classifier(
      k = 40,
      use_approx_nn = TRUE,
      distance_metric = "euclidean",
      exp_doub_rate = 0.1,
      stdev_doub_rate = 0.03,
      get_neighbor_parents = FALSE
    )

------------------------------------------------------------------------

### Method `plot_histogram()`

#### Usage

    ScrubletR$plot_histogram(
      scale_hist_obs = "log",
      scale_hist_sim = "linear",
      fig_size = c(8, 3)
    )

------------------------------------------------------------------------

### Method `set_embedding()`

#### Usage

    ScrubletR$set_embedding(embedding_name, coordinates)

------------------------------------------------------------------------

### Method `plot_embedding()`

#### Usage

    ScrubletR$plot_embedding(
      embedding_name,
      score = "raw",
      marker_size = 5,
      order_points = FALSE,
      fig_size = c(8, 4),
      color_map = NULL
    )

------------------------------------------------------------------------

### Method `pipeline_normalize()`

#### Usage

    ScrubletR$pipeline_normalize(postnorm_total = NULL)

------------------------------------------------------------------------

### Method `pipeline_get_gene_filter()`

#### Usage

    ScrubletR$pipeline_get_gene_filter(
      min_counts = 3,
      min_cells = 3,
      min_gene_variability_pctl = 85,
      plot = T
    )

------------------------------------------------------------------------

### Method `pipeline_apply_gene_filter()`

#### Usage

    ScrubletR$pipeline_apply_gene_filter()

------------------------------------------------------------------------

### Method `pipeline_mean_center()`

#### Usage

    ScrubletR$pipeline_mean_center()

------------------------------------------------------------------------

### Method `pipeline_normalize_variance()`

#### Usage

    ScrubletR$pipeline_normalize_variance()

------------------------------------------------------------------------

### Method `pipeline_zscore()`

#### Usage

    ScrubletR$pipeline_zscore()

------------------------------------------------------------------------

### Method `pipeline_log_transform()`

#### Usage

    ScrubletR$pipeline_log_transform(pseudocount = 1)

------------------------------------------------------------------------

### Method `pipeline_truncated_svd()`

#### Usage

    ScrubletR$pipeline_truncated_svd(
      n_prin_comps = 30,
      random_state = 0,
      algorithm = "arpack"
    )

------------------------------------------------------------------------

### Method `pipeline_pca()`

#### Usage

    ScrubletR$pipeline_pca(
      n_prin_comps = 50,
      random_state = 0,
      svd_solver = "arpack"
    )

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    ScrubletR$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
