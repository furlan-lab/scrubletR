# Initialize Scrublet object with counts matrix and doublet prediction parameters
#
#         Parameters
#         ----------
#         counts_matrix : scipy sparse matrix or ndarray, shape (n_cells, n_genes)
#             Matrix containing raw (unnormalized) UMI-based transcript counts.
#             Converted into a scipy.sparse.csc_matrix.
#
#         total_counts : ndarray, shape (n_cells,), optional (default: None)
#             Array of total UMI counts per cell. If `None`, this is calculated
#             as the row sums of `counts_matrix`.
#
#         sim_doublet_ratio : float, optional (default: 2.0)
#             Number of doublets to simulate relative to the number of observed
#             transcriptomes.
#
#         n_neighbors : int, optional (default: None)
#             Number of neighbors used to construct the KNN graph of observed
#             transcriptomes and simulated doublets. If `None`, this is
#             set to round(0.5 * sqrt(n_cells))
#
#         expected_doublet_rate : float, optional (default: 0.1)
#             The estimated doublet rate for the experiment.
#
#         stdev_doublet_rate : float, optional (default: 0.02)
#             Uncertainty in the expected doublet rate.
#
#         random_state : int, optional (default: 0)
#             Random state for doublet simulation, approximate
#             nearest neighbor search, and PCA/TruncatedSVD.
#
#         Attributes
#         ----------
#         predicted_doublets_ : ndarray, shape (n_cells,)
#             Boolean mask of predicted doublets in the observed
#             transcriptomes.
#
#         doublet_scores_obs_ : ndarray, shape (n_cells,)
#             Doublet scores for observed transcriptomes.
#
#         doublet_scores_sim_ : ndarray, shape (n_doublets,)
#             Doublet scores for simulated doublets.
#
#         doublet_errors_obs_ : ndarray, shape (n_cells,)
#             Standard error in the doublet scores for observed
#             transcriptomes.
#
#         doublet_errors_sim_ : ndarray, shape (n_doublets,)
#             Standard error in the doublet scores for simulated
#             doublets.
#
#         threshold_: float
#             Doublet score threshold for calling a transcriptome
#             a doublet.
#
#         z_scores_ : ndarray, shape (n_cells,)
#             Z-score conveying confidence in doublet calls.
#             Z = `(doublet_score_obs_ - threhsold_) / doublet_errors_obs_`
#
#         detected_doublet_rate_: float
#             Fraction of observed transcriptomes that have been called
#             doublets.
#
#         detectable_doublet_fraction_: float
#             Estimated fraction of doublets that are detectable, i.e.,
#             fraction of simulated doublets with doublet scores above
#             `threshold_`
#
#         overall_doublet_rate_: float
#             Estimated overall doublet rate,
#             `detected_doublet_rate_ / detectable_doublet_fraction_`.
#             Should agree (roughly) with `expected_doublet_rate`.
#
#         manifold_obs_: ndarray, shape (n_cells, n_features)
#             The single-cell "manifold" coordinates (e.g., PCA coordinates)
#             for observed transcriptomes. Nearest neighbors are found using
#             the union of `manifold_obs_` and `manifold_sim_` (see below).
#
#         manifold_sim_: ndarray, shape (n_doublets, n_features)
#             The single-cell "manifold" coordinates (e.g., PCA coordinates)
#             for simulated doublets. Nearest neighbors are found using
#             the union of `manifold_obs_` (see above) and `manifold_sim_`.
#
#         doublet_parents : ndarray, shape (n_doublets, 2)
#             Indices of the observed transcriptomes used to generate the
#             simulated doublets.
#
#         doublet_neighbor_parents_ : list, length n_cells
#             A list of arrays of the indices of the doublet neighbors of
#             each observed transcriptome (the ith entry is an array of
#             the doublet neighbors of transcriptome i).


#' Scrublet R6 Class
#'
#' @importFrom R6 R6Class
#' @import Matrix
#' @export

Scrublet <- R6::R6Class("Scrublet",
                    public = list(
                      E_obs = NULL,
                      E_sim = NULL,
                      E_obs_norm = NULL,
                      E_sim_norm = NULL,
                      gene_filter = NULL,
                      embeddings = NULL,
                      total_counts_obs = NULL,
                      total_counts_sim = NULL,
                      sim_doublet_ratio = NULL,
                      n_neighbors = NULL,
                      expected_doublet_rate = NULL,
                      stdev_doublet_rate = NULL,
                      random_state = NULL,
                      doublet_parents = NULL,
                      manifold_obs_ = NULL,
                      manifold_sim_ = NULL,
                      doublet_scores_obs_ = NULL,
                      doublet_scores_sim_ = NULL,
                      doublet_errors_obs_ = NULL,
                      doublet_errors_sim_ = NULL,
                      doublet_neighbor_parents_ = NULL,
                      predicted_doublets = NULL,
                      z_scores_ = NULL,
                      threshold_ = NULL,
                      detected_doublet_rate_ = NULL,
                      detectable_doublet_fraction_ = NULL,
                      overall_doublet_rate_ = NULL,
                      get_dims = function(){
                        message(paste0("E_obs: ", dim(self$E_obs)[1], " rows by ", dim(self$E_obs)[2], " columns"))
                        message(paste0("E_sim: ", dim(self$E_sim)[1], " rows by ", dim(self$E_sim)[2], " columns"))
                        message(paste0("E_obs_norm: ", dim(self$E_obs_norm)[1], " rows by ", dim(self$E_obs_norm)[2], " columns"))
                        message(paste0("E_sim_norm: ", dim(self$E_sim_norm)[1], " rows by ", dim(self$E_sim_norm)[2], " columns"))
                      },
                      initialize = function(
                                            counts_matrix,
                                            total_counts = NULL,
                                            sim_doublet_ratio = 2.0,
                                            n_neighbors = NULL,
                                            expected_doublet_rate = 0.1,
                                            stdev_doublet_rate = 0.02,
                                            random_state = 0
                      ) {
                        # Convert counts_matrix to a sparse matrix if not already
                        # if (!is.sparse(counts_matrix)) {
                        #   counts_matrix <- Matrix::Matrix(counts_matrix, sparse = TRUE)
                        # } else if (!is(CscMatrix, counts_matrix)) {
                        #   counts_matrix <- as(counts_matrix, "CsparseMatrix")
                        # }

                        # Initialize counts matrices
                        self$E_obs <- counts_matrix
                        self$E_sim <- NULL
                        self$E_obs_norm <- NULL
                        self$E_sim_norm <- NULL

                        # Initialize total_counts_obs
                        if (is.null(total_counts)) {
                          self$total_counts_obs <- Matrix::rowSums(self$E_obs)
                        } else {
                          self$total_counts_obs <- total_counts
                        }

                        # Initialize gene filter
                        self$gene_filter <- 1:ncol(self$E_obs)

                        # Initialize embeddings
                        self$embeddings <- list()

                        # Set other attributes
                        self$sim_doublet_ratio = sim_doublet_ratio
                        self$n_neighbors = n_neighbors
                        self$expected_doublet_rate = expected_doublet_rate
                        self$stdev_doublet_rate = stdev_doublet_rate
                        self$random_state = random_state
                        # Set n_neighbors if NULL
                        if (is.null(self$n_neighbors)) {
                          self$n_neighbors <- round(0.5 * sqrt(nrow(self$E_obs)))
                        }
                      },

    scrub_doublets = function(
                              synthetic_doublet_umi_subsampling = 1.0,
                              use_approx_neighbors = TRUE,
                              distance_metric = 'euclidean',
                              get_doublet_neighbor_parents = FALSE,
                              min_counts = 3,
                              min_cells = 3,
                              min_gene_variability_pctl = 85,
                              log_transform = FALSE,
                              mean_center = T,
                              normalize_variance = T,
                              n_prin_comps = 30,
                              svd_solver = 'arpack',
                              verbose = TRUE
                              ) {
      t0 <- Sys.time()

      self$E_sim <- NULL
      self$E_obs_norm <- NULL
      self$E_sim_norm <- NULL
      self$gene_filter <- 1:ncol(self$E_obs)

      cat_optional('Preprocessing...', verbose)
      self$pipeline_normalize()
      self$pipeline_get_gene_filter(min_counts=min_counts, min_cells=min_cells, min_gene_variability_pctl=min_gene_variability_pctl)
      self$pipeline_apply_gene_filter()

      cat_optional('Simulating doublets...', verbose)
      self$simulate_doublets(sim_doublet_ratio=self$sim_doublet_ratio, synthetic_doublet_umi_subsampling=synthetic_doublet_umi_subsampling)
      self$pipeline_normalize(postnorm_total=1e6)

      if (log_transform) {
        self$pipeline_log_transform()
      }

      if (mean_center && normalize_variance) {
        self$pipeline_zscore()
      } else if (mean_center) {
        self$pipeline_mean_center()
      } else if (normalize_variance) {
        self$pipeline_normalize_variance()
      }

      if (mean_center) {
        cat_optional('Embedding transcriptomes using PCA...', verbose)
        self$pipeline_pca( n_prin_comps=n_prin_comps, random_state=self$random_state, svd_solver=svd_solver)
      } else {
        cat_optional('Embedding transcriptomes using Truncated SVD...', verbose)
        self$pipeline_truncated_svd( n_prin_comps=n_prin_comps, random_state=self$random_state, algorithm=svd_solver)
      }

      cat_optional('Calculating doublet scores...', verbose)
      self$calculate_doublet_scores(use_approx_neighbors=use_approx_neighbors,
                               distance_metric=distance_metric,
                               get_doublet_neighbor_parents=get_doublet_neighbor_parents)
      self$call_doublets(verbose=verbose)

      t1 <- Sys.time()
      cat_optional(sprintf('Elapsed time: %.1f seconds', as.numeric(difftime(t1, t0, units = "secs"))), verbose)
      return(list(doublet_scores = self$doublet_scores_obs_, predicted_doublets = self$predicted_doublets))
    },

    simulate_doublets = function(
                                  sim_doublet_ratio = NULL,
                                  synthetic_doublet_umi_subsampling = 1.0
                                  ) {
      # Simulate doublets by adding the counts of random observed transcriptome pairs.

      # Arguments
      # sim_doublet_ratio: float, optional (default: NULL)
      #   Number of doublets to simulate relative to the number of observed transcriptomes.
      #   If NULL, self$sim_doublet_ratio is used.
      # synthetic_doublet_umi_subsampling: float, optional (default: 1.0)
      #   Rate for sampling UMIs when creating synthetic doublets.
      #   If 1.0, each doublet is created by simply adding the UMIs from two randomly
      #   sampled observed transcriptomes. For values less than 1, the UMI counts are
      #   added and then randomly sampled at the specified rate.

      # Sets
      # doublet_parents

      if (is.null(sim_doublet_ratio)) {
        sim_doublet_ratio = self$sim_doublet_ratio
      } else {
        self$sim_doublet_ratio = sim_doublet_ratio
      }

      n_obs <- nrow(self$E_obs)
      n_sim <- as.integer(n_obs * sim_doublet_ratio)

      set.seed(self$random_state)
      pair_ix <- matrix(sample(1:n_obs, size = n_sim * 2, replace = TRUE), ncol = 2)

      E1 <- self$E_obs[pair_ix[, 1], , drop = FALSE]
      E2 <- self$E_obs[pair_ix[, 2], , drop = FALSE]
      tots1 <- self$total_counts_obs[pair_ix[, 1]]
      tots2 <- self$total_counts_obs[pair_ix[, 2]]

      if (synthetic_doublet_umi_subsampling < 1) {
        result <- subsample_counts(E1 + E2, synthetic_doublet_umi_subsampling, tots1 + tots2, random_seed = self$random_state)
        self$E_sim <- result$counts
        self$total_counts_sim <- result$total_counts
      } else {
        self$E_sim <- E1 + E2
        self$total_counts_sim <- tots1 + tots2
      }

      self$doublet_parents <- pair_ix
      return(NULL)
    },

    set_manifold = function(manifold_obs, manifold_sim) {
      self$manifold_obs_ = manifold_obs
      self$manifold_sim_ = manifold_sim
    },

    calculate_doublet_scores = function(
                                        use_approx_neighbors = TRUE,
                                        distance_metric = 'euclidean',
                                        get_doublet_neighbor_parents = FALSE
                                        ) {
      self$nearest_neighbor_classifier(
        k=self$n_neighbors,
        exp_doub_rate=self$expected_doublet_rate,
        stdev_doub_rate=self$stdev_doublet_rate,
        use_approx_nn=use_approx_neighbors,
        distance_metric=distance_metric,
        get_neighbor_parents=get_doublet_neighbor_parents
      )
    },

    call_doublets = function(threshold = NULL, verbose = TRUE) {
      if (is.null(threshold)) {
        # Automatic threshold detection
        tryCatch({
          threshold <- threshold_minimum(self$doublet_scores_sim_)
          cat_optional(paste0("Automatically set threshold at doublet score =", threshold), verbose)
        }, error = function(e) {
          self$predicted_doublets <- NULL
          cat_optional(paste0("Warning: failed to automatically identify doublet score threshold. Run `call_doublets` with a user-specified threshold."), verbose)
          return(self$predicted_doublets)
        })
      }
      Ld_obs <- self$doublet_scores_obs_
      Ld_sim <- self$doublet_scores_sim_
      se_obs <- self$doublet_errors_obs_
      Z <- (Ld_obs - threshold) / se_obs
      self$predicted_doublets <- Ld_obs > threshold
      self$z_scores_ <- Z
      self$threshold_ <- threshold
      self$detected_doublet_rate_ <- sum(Ld_obs > threshold) / length(Ld_obs)
      self$detectable_doublet_fraction_ <- sum(Ld_sim > threshold) / length(Ld_sim)
      self$overall_doublet_rate_ <- self$detected_doublet_rate_ / self$detectable_doublet_fraction_
      cat_optional(paste0('Detected doublet rate =', 100 * self$detected_doublet_rate_, "%"), verbose)
      cat_optional(paste0('Estimated detectable doublet fraction =', 100 * self$detectable_doublet_fraction_, "%"), verbose)
      cat_optional(paste0('Overall doublet rate:'), verbose)
      cat_optional(paste0('\tExpected   =', 100 * self$expected_doublet_rate, "%"), verbose)
      cat_optional(paste0('\tEstimated  =', 100 * self$overall_doublet_rate_, "%"), verbose)
    },

    #' Nearest Neighbor Classifier
    #'
    #' This method performs nearest neighbor classification to calculate doublet scores.
    #'
    #' @param k Number of nearest neighbors.
    #' @param use_approx_nn Whether to use approximate nearest neighbors.
    #' @param distance_metric Distance metric for finding neighbors.
    #' @param exp_doub_rate Expected doublet rate.
    #' @param stdev_doub_rate Standard deviation of doublet rate.
    #' @param get_neighbor_parents Whether to get parents of doublet neighbors.
    #'
    #' @return NULL (modifies the R6 object in place)
    nearest_neighbor_classifier = function(
                                          k = 40,
                                          use_approx_nn = TRUE,
                                          distance_metric = 'euclidean',
                                          exp_doub_rate = 0.1,
                                          stdev_doub_rate = 0.03,
                                          get_neighbor_parents = FALSE) {
      manifold <- rbind(self$manifold_obs_, self$manifold_sim_)
      doub_labels <- c(rep(0, nrow(self$manifold_obs_)), rep(1, nrow(self$manifold_sim_)))

      n_obs <- sum(doub_labels == 0)
      n_sim <- sum(doub_labels == 1)

      # Adjust k (number of nearest neighbors) based on the ratio of simulated to observed cells
      k_adj <- round(k * (1 + n_sim / n_obs))

      # Find k_adj nearest neighbors
      neighbors <- get_knn_graph(manifold, k = k_adj, dist_metric = distance_metric, approx = use_approx_nn, return_edges = FALSE, random_seed = self$random_state)
      #dim(neighbors)
      # Calculate doublet score based on the ratio of simulated cell neighbors vs. observed cell neighbors
      # doub_neigh_mask <- doub_labels[neighbors] == 1
      # n_sim_neigh <- rowSums(doub_neigh_mask)
      # n_obs_neigh <- ncol(doub_neigh_mask) - n_sim_neigh

      doub_neigh_mask <- matrix(doub_labels[neighbors] == 1, nrow = n_obs + n_sim)

      n_sim_neigh <- apply(doub_neigh_mask, 1, sum)
      n_obs_neigh <- ncol(doub_neigh_mask) - n_sim_neigh

      rho <- exp_doub_rate
      r <- n_sim / n_obs
      nd <- as.numeric(n_sim_neigh)
      ns <- as.numeric(n_obs_neigh)
      N <- as.numeric(k_adj)

      # Bayesian
      q <- (nd + 1) / (N + 2)
      Ld <- q * rho / r / (1 - rho - q * (1 - rho - rho / r))

      se_q <- sqrt(q * (1 - q) / (N + 3))
      se_rho <- stdev_doub_rate

      se_Ld <- q * rho / r / (1 - rho - q * (1 - rho - rho / r))^2 * sqrt((se_q / q * (1 - rho))^2 + (se_rho / rho * (1 - q))^2)

      self$doublet_scores_obs_ <- Ld[doub_labels == 0]
      self$doublet_scores_sim_ <- Ld[doub_labels == 1]
      self$doublet_errors_obs_ <- se_Ld[doub_labels == 0]
      self$doublet_errors_sim_ <- se_Ld[doub_labels == 1]

      # get parents of doublet neighbors, if requested
      neighbor_parents <- NULL
      if (get_neighbor_parents) {
        parent_cells <- self$doublet_parents_
        neighbors <- neighbors - n_obs
        neighbor_parents <- list()

        for (iCell in seq_len(n_obs)) {
          this_doub_neigh <- neighbors[iCell, neighbors[iCell, ] > -1]

          if (length(this_doub_neigh) > 0) {
            this_doub_neigh_parents <- unique(parent_cells[this_doub_neigh, , drop = FALSE])
            neighbor_parents[[iCell]] <- this_doub_neigh_parents
          } else {
            neighbor_parents[[iCell]] <- integer(0)
          }
        }

        self$doublet_neighbor_parents_ <- neighbor_parents
      }
    },


    plot_histogram = function(
                              scale_hist_obs = 'log',
                              scale_hist_sim = 'linear',
                              fig_size = c(8, 3)
                              ) {
      # Plot histogram code here
    },

    set_embedding = function(embedding_name, coordinates) {
      # Set embedding code here
    },

    plot_embedding = function(
                            embedding_name,
                            score = 'raw',
                            marker_size = 5,
                            order_points = FALSE,
                            fig_size = c(8, 4),
                            color_map = NULL
                            ) {
      # Plot embedding code here
    },

    #' Total Counts Normalization
    #'
    #' This function performs total counts normalization on the observed and simulated expression matrices.
    #'
    #' @param postnorm_total Target total count after normalization. If not provided, the mean of observed total counts will be used.
    #'
    #' @return NULL (modifies the R6 object in place)
    #'
    #' @details
    #' Total counts normalization is a step in the preprocessing of single-cell RNA-seq data. It aims to scale the total counts of each cell so that they have a consistent value, making them comparable across cells. This is useful for downstream analyses, such as dimensionality reduction and clustering.
    #'
    #' The function updates the normalized observed expression matrix (`E_obs_norm`) based on the specified target total count. If a simulated expression matrix (`E_sim`) is present in the R6 object, it also normalizes the simulated matrix (`E_sim_norm`).
    #'
    #' @seealso
    #' \code{\link{tot_counts_norm}}

    pipeline_normalize = function(postnorm_total = NULL) {
        # Total counts normalization
        if (is.null(postnorm_total)) {
          postnorm_total = mean(self$total_counts_obs)
        }

        self$E_obs_norm <- tot_counts_norm(self$E_obs, target_total = postnorm_total, total_counts = self$total_counts_obs)

        if (!is.null(self$E_sim)) {
          self$E_sim_norm <- tot_counts_norm(self$E_sim, target_total = postnorm_total, total_counts = self$total_counts_sim)
        }
      },
    #' Identify highly variable genes expressed above a minimum level
    #' @param min_counts Minimum counts for a gene to be considered highly variable
    #' @param min_cells Minimum cells for a gene to be considered highly variable
    #' @param min_vscore_pctl Minimum percentile of variance score for a gene to be considered highly variable
    #' @return NULL (modifies the R6 object in place)
    pipeline_get_gene_filter = function(
                                        min_counts = 3,
                                        min_cells = 3,
                                        min_gene_variability_pctl = 85
                                        ) {
      # counts <- Matrix::rowSums(self$E_obs_norm)
      # cells <- Matrix::colSums(self$E_obs_norm > 0)
      # vscores <- get_vscores(self$E_obs_norm)
      #
      # highly_variable <- counts >= min_counts & cells >= min_cells & vscores >= quantile(vscores, min_vscore_pctl / 100)
      # return(highly_variable)
      self$gene_filter = filter_genes(self$E_obs_norm,
                                       min_counts=min_counts,
                                       min_cells=min_cells,
                                       min_vscore_pctl=min_gene_variability_pctl)
    },



    #' Apply gene filter to expression matrices
    #' @return NULL (modifies the R6 object in place)
    pipeline_apply_gene_filter = function() {
      self$E_obs <- self$E_obs[, self$gene_filter, drop = FALSE]
      self$E_obs_norm <- self$E_obs_norm[, self$gene_filter, drop = FALSE]
      self$E_sim <- self$E_sim[, self$gene_filter, drop = FALSE]
      self$E_sim_norm <- self$E_sim_norm[, self$gene_filter, drop = FALSE]
    },

    #' Mean center expression matrix
    #' @return NULL (modifies the R6 object in place)
    pipeline_mean_center = function() {
      gene_means <- colMeans(self$E_obs_norm, na.rm = TRUE)
      self$E_obs_norm <- sweep(self$E_obs_norm, 2, gene_means, "-")
      if (!is.null(self$E_sim_norm)) {
        self$E_sim_norm <- sweep(self$E_sim_norm, 2, gene_means, "-")
      }
    },

    #' Variance Normalization of Expression Matrices
    #'
    #' This function performs variance normalization on the observed and simulated expression matrices.
    #'
    #' @return NULL (modifies the R6 object in place)
    #'
    #' @details
    #' Variance normalization is a step that scales the expression values of each gene across cells based on their standard deviation. It helps to make gene expression distributions comparable between genes.
    #'
    #' The function calculates gene standard deviations from the observed expression matrix (`E_obs_norm`). It then scales the observed and simulated matrices based on these standard deviations.


    pipeline_normalize_variance = function() {
      gene_stdevs <- sqrt(sparse_var(self$E_obs_norm))
      self$E_obs_norm <- t(t(self$E_obs_norm) / gene_stdevs)
      if (!is.null(self$E_sim_norm)) {
        self$E_sim_norm <- t(t(self$E_sim_norm) / gene_stdevs)
      }
    },

    #' Z-Score Normalization of Expression Matrices
    #'
    #' This function performs z-score normalization on the observed and simulated expression matrices.
    #'
    #' @return NULL (modifies the R6 object in place)
    #' @importFrom Matrix rowMeans
    #' @details
    #' Z-score normalization standardizes the expression values of each gene across cells by subtracting the mean and dividing by the standard deviation. It helps to make gene expression distributions comparable between genes.
    #'
    #' The function calculates gene means and standard deviations from the observed expression matrix (`E_obs_norm`). It then applies z-score normalization to the observed and simulated matrices.

    pipeline_zscore = function() {
      gene_means = colMeans(self$E_obs_norm) #######changed to column??? SCOTT
      gene_stdevs = sqrt(sparse_var(self$E_obs_norm, axis = 2))
      self$E_obs_norm = sparse_zscore(self$E_obs_norm, gene_means, gene_stdevs)
      if (!is.null(self$E_sim_norm)) {
        self$E_sim_norm <- sparse_zscore(self$E_sim_norm, gene_means, gene_stdevs)
      }
    },

    #' Log Transform Expression Matrices
    #'
    #' This method performs log transformation on the observed and simulated expression matrices.
    #'
    #' @param pseudocount The pseudocount added to expression values before taking the logarithm.
    #'
    #' @return NULL (modifies the R6 object in place)

    pipeline_log_transform = function(pseudocount = 1) {
      self$E_obs_norm <- log_normalize(self$E_obs_norm, pseudocount)

      if (!is.null(self$E_sim_norm)) {
        self$E_sim_norm <- log_normalize(self$E_sim_norm, pseudocount)
      }

      invisible(NULL)
    },

    #' Truncated Singular Value Decomposition (SVD) for Expression Matrices
    #'
    #' This method performs truncated SVD on the observed expression matrix and sets the resulting manifold.
    #'
    #' @param n_prin_comps Number of principal components to keep.
    #' @param random_state Seed for reproducibility.
    #' @param algorithm Algorithm to use for SVD.
    #'
    #' @return NULL (modifies the R6 object in place)
    pipeline_truncated_svd = function(n_prin_comps = 30, random_state = 0, algorithm = 'arpack') {
      svd <- svd(self$E_obs_norm, nu = n_prin_comps, nv = 0)
      ####### not sure about this next line SCOTT
      #self$set_manifold(predict(svd, newdata = self$E_obs_norm), predict(svd, newdata = self$E_sim_norm))
    },

    #' Principal Component Analysis (PCA) for Expression Matrices
    #'
    #' This method performs PCA on the observed expression matrix and sets the resulting manifold.
    #'
    #' @param n_prin_comps Number of principal components to keep.
    #' @param random_state Seed for reproducibility.
    #' @param svd_solver Algorithm to use for SVD.
    #' @importFrom irlba prcomp_irlba
    #' @return NULL (modifies the R6 object in place)

    pipeline_pca = function(n_prin_comps = 50, random_state = 0, svd_solver = 'arpack') {
      X_obs <- as.matrix(self$E_obs_norm)

      if (!is.null(self$E_sim_norm)) {
        X_sim <- as.matrix(self$E_sim_norm)
      } else {
        X_sim <- NULL
      }
      pca <- prcomp_irlba(X_obs, n = n_prin_comps, center = TRUE, scale. = FALSE)
      # pca <- prcomp(X_obs, center = TRUE, scale. = FALSE)
      self$set_manifold(pca$x[, 1:n_prin_comps], pca$rotation[, 1:n_prin_comps])

    }
                    )
)
