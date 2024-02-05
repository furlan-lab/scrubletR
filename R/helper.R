#' Total Counts Normalization
#' @param E Counts matrix to be normalized
#' @param total_counts Vector of total counts per cell (if NULL, computed from E)
#' @param exclude_dominant_frac Exclude overly abundant genes fraction (default: 1)
#' @param included Indices of genes to include in normalization (if empty, all genes are considered)
#' @param target_total Target total for normalization (if NULL, the mean of total counts is used)
#' @return List containing normalized counts matrix (Enorm), average total counts, and vector of included genes
tot_counts_norm <- function(E, total_counts = NULL, exclude_dominant_frac = 1, included = integer(0), target_total = NULL) {
  ncell <- nrow(E)
  # E <- as(E, "Matrix")

  if (is.null(total_counts)) {
    if (length(included) == 0) {
      if (exclude_dominant_frac == 1) {
        tots_use <- Matrix::colSums(E)
      } else {
        tots <- Matrix::colSums(E)
        wtmp <- Matrix::Matrix(0, ncell, ncell)
        Matrix::diag(wtmp) <- 1 / tots
        included <- which(!(Matrix::colSums(wtmp %*% E > exclude_dominant_frac) > 0))
        tots_use <- Matrix::colSums(E[, included])
        cat('Excluded', sum(!included), 'genes from normalization\n')
      }
    } else {
      tots_use <- Matrix::colSums(E[, included])
    }
  } else {
    tots_use <- total_counts
  }

  if (is.null(target_total)) {
    target_total <- mean(tots_use)
  }

  w <- Matrix::Matrix(0, ncell, ncell)
  Matrix::diag(w) <- target_total / tots_use
  Enorm <- w %*% E

  return(Enorm)
}

#' Custom messaging
#' @param message message to be printed if verbose
#' @param verbose verbosity (bool)
#' @return NULL
#' @export
cat_optional <- function(message, verbose) {
  if (verbose) {
    message(paste0(message))
  }
}


#' Get k-Nearest-Neighbor Graph
#'
#' Build k-nearest-neighbor graph and return edge list and nearest neighbor matrix.
#'
#' @param X Data matrix.
#' @param k Number of nearest neighbors.
#' @param dist_metric Distance metric for finding neighbors.
#' @param approx Whether to use approximate nearest neighbor search.
#' @param return_edges Whether to return edge list.
#' @param random_seed Random seed for reproducibility.
#' @import RcppAnnoy
#' @importFrom methods new
#' @return If \code{return_edges} is \code{TRUE}, a list containing edge list and nearest neighbor matrix; otherwise, the nearest neighbor matrix.
#' @export
get_knn_graph <- function(X, k = 5, dist_metric = c('euclidean', 'angular', 'manhattan', 'hamming'), approx = FALSE, return_edges = TRUE, random_seed = 0) {
  t0 <- Sys.time()
  dist_metric <- match.arg(dist_metric)

  if (approx) {
    if (!requireNamespace("RcppAnnoy", quietly = TRUE)) {
      approx <- FALSE
      cat('Could not find library "annoy" for approx. nearest neighbor search\n')
    }
  }

  if (approx) {
    npc <- ncol(X)
    ncell <- nrow(X)
    switch (dist_metric,
      euclidean = {annoy_class = AnnoyEuclidean},
      angular = {annoy_class = AnnoyAngular},
      manhattan = {annoy_class = AnnoyManhattan},
      hamming = {annoy_class = AnnoyHamming},
    )
    a <- new(annoy_class, npc)
    a$setSeed(random_seed)

    for (i in seq_len(ncell)) {
      a$addItem(i - 1, as.numeric(X[i, , drop = FALSE]))
    }
    a$build(10)  # 10 trees

    knn <- lapply(seq_len(ncell), function(iCell) {
      a$getNNsByItem(iCell - 1, k + 1)[-1]
    })
    knn <- do.call(rbind, knn) + 1
  } else {
    # if (dist_metric == 'cosine') {
    #   nbrs <- get.knn(X, k = k, method = 'cosine')
    # } else {
    #   nbrs <- get.knn(X, k = k, method = 'euclidean')
    # }
    # knn <- nbrs$nn.index + 1
    stop("not implemented yet")
  }

  if (return_edges) {
    links <- lapply(seq_len(nrow(knn)), function(i) {
      sort(c(i, knn[i, ]))
    })
    links <- unique(links)

    t_elapse <- difftime(Sys.time(), t0, units = "secs")
    cat(sprintf('kNN graph built in %.3f sec\n', t_elapse))

    return(list(links = links, knn = knn))
  }
  return(knn)
}

#' Build Adjacency Matrix
#'
#' Build an adjacency matrix from the provided edge list.
#'
#' @param edges Edge list.
#' @param n_nodes Number of nodes.
#'
#' @return Adjacency matrix.
#' @importFrom Matrix sparseMatrix
#' @export
build_adj_mat <- function(edges, n_nodes) {
  A <- sparseMatrix(i = sapply(edges, function(e) e[1] - 1), j = sapply(edges, function(e) e[2] - 1), x = rep(1, length(edges)), dims = c(n_nodes, n_nodes))
  return(A)
}



#' Calculate v-scores and related statistics for genes in the input counts matrix
#'
#' This function calculates v-scores, coefficient of variation (CV), and other statistics
#' for genes in the input counts matrix. The v-score is an above-Poisson noise statistic
#' that helps assess the variability of gene expression levels.
#'
#' @param E A counts matrix where rows represent cells and columns represent genes.
#' @param min_mean Minimum mean expression value for genes to be considered.
#' @param nBins Number of bins for calculating running quantiles.
#' @param fit_percentile Percentile used for fitting the running quantile.
#' @param error_wt Weight for the error function during optimization.
#'
#' @return A list containing the following components:
#'   \itemize{
#'     \item \code{v_scores}: V-scores for each gene.
#'     \item \code{CV_eff}: Coefficient of variation (CV) for effective counts.
#'     \item \code{CV_input}: Coefficient of variation (CV) for input counts.
#'     \item \code{gene_ix}: Indices of genes considered in the analysis.
#'     \item \code{mu_gene}: Mean expression values for selected genes.
#'     \item \code{FF_gene}: Fano factor (variance to mean ratio) for selected genes.
#'     \item \code{a}: Parameter 'a' obtained during optimization.
#'     \item \code{b}: Parameter 'b' obtained during optimization.
#'   }
#'
#' @seealso \code{\link{runningquantile}}, \code{\link{optimize}}
#'
#' @importFrom graphics hist
#' @importFrom stats optimise
#' @export
get_vscores <- function(E, min_mean = 0, nBins = 50, fit_percentile = 0.1, error_wt = 1) {
  ncell <- nrow(E)
  mu_gene <- colMeans(E)
  gene_ix <- which(mu_gene > min_mean)
  mu_gene <- mu_gene[gene_ix]

  tmp <- E[, gene_ix]
  tmp <- tmp^2
  var_gene <- colMeans(tmp) - mu_gene^2
  rm(tmp)
  FF_gene <- var_gene / mu_gene

  data_x <- log(mu_gene)
  data_y <- log(FF_gene / mu_gene)


  result <- runningquantile(data_x, data_y, fit_percentile, nBins)
  x <- result$xOut[!is.na(result$yOut)]
  y <- result$yOut[!is.na(result$yOut)]

  gLog <- function(input1, input2, input3) {
    log(input2 * exp(-input1) + input3)
  }

  hist_data <- hist(log(FF_gene[mu_gene > 0]), breaks = 200, plot = F)
  b <- hist_data$breaks[-length(hist_data$breaks)] + diff(hist_data$breaks) / 2
  max_ix <- which.max(hist_data$counts)
  c <- max(exp(b[max_ix]), 1)

  errFun <- function(b2) {
    sum(abs(gLog(x, c, b2) - y)^error_wt)
  }


  b0 <- 0.1
  ###CHANGED THIS SCOTT
  b<- optimise(errFun, c(b0, 1e30*b0))$minimum
  a <- c / (1 + b) - 1

  v_scores <- FF_gene / ((1 + a) * (1 + b) + b * mu_gene)
  CV_eff <- sqrt((1 + a) * (1 + b) - 1)
  CV_input <- sqrt(b)

  result <- list(v_scores = v_scores,
                 CV_eff = CV_eff,
                 CV_input = CV_input,
                 gene_ix = gene_ix,
                 mu_gene = mu_gene,
                 FF_gene = FF_gene,
                 a = a,
                 b = b)

  return(result)
}


#' Calculate running quantiles for given data points
#'
#' This function calculates running quantiles for given x and y data points.
#'
#' @param x Numeric vector representing the x-axis data points.
#' @param y Numeric vector representing the y-axis data points.
#' @param p Percentile value for calculating the running quantiles.
#' @param nBins Number of bins for the running quantiles.
#'
#' @return A list containing the following components:
#'   \itemize{
#'     \item \code{xOut}: x-axis values for the running quantiles.
#'     \item \code{yOut}: y-axis values representing the running quantiles.
#'   }
#' @importFrom stats quantile
#' @export
runningquantile <- function(x, y, p, nBins) {
  ind <- order(x)
  x <- x[ind]
  y <- y[ind]

  dx <- (x[length(x)] - x[1]) / nBins
  xOut <- seq(x[1] + dx/2, x[length(x)] - dx/2, length.out = nBins)

  yOut <- numeric(length = length(xOut))

  for (i in seq_along(xOut)) {
    ind <- which(x >= xOut[i] - dx/2 & x < xOut[i] + dx/2)
    if (length(ind) > 0) {
      yOut[i] <- quantile(y[ind], p)
    } else {
      if (i > 1) {
        yOut[i] <- yOut[i - 1]
      } else {
        yOut[i] <- NaN
      }
    }
  }

  result <- list(xOut = xOut, yOut = yOut)
  return(result)
}

#' Filter genes by expression level and variability
#'
#' This function filters genes based on expression level and variability using v-scores.
#'
#' @param E A counts matrix where rows represent cells and columns represent genes.
#' @param base_ix Indices of cells to be used for v-score calculation (default: all cells).
#' @param min_vscore_pctl Minimum percentile threshold for v-scores.
#' @param min_counts Minimum expression counts required for a gene to be considered.
#' @param min_cells Minimum number of cells expressing a gene for it to be considered.
#' @param show_vscore_plot Logical indicating whether to show a v-score plot.
#' @param sample_name Character string specifying the name of the sample for the plot title.
#'
#' @return A numeric vector containing the indices of filtered genes.
#' @seealso \code{\link{get_vscores}}
#' @importFrom stats quantile
#' @importFrom grDevices rgb
#' @importFrom graphics points
#' @importFrom graphics lines
#' @import ggplot2

filter_genes <- function(E, base_ix = NULL, min_vscore_pctl = 85, min_counts = 3, min_cells = 3,
                         plot = FALSE, sample_name = '') {
  if (is.null(base_ix)) {
    base_ix <- seq_len(nrow(E))
  }

  vscores_result <- get_vscores(E[base_ix, , drop = FALSE])
  Vscores <- vscores_result$v_scores
  gene_ix <- vscores_result$gene_ix
  mu_gene <- vscores_result$mu_gene
  FF_gene <- vscores_result$FF_gene
  a <- vscores_result$a
  b <- vscores_result$b

  ix2 <- Vscores > 0
  Vscores <- Vscores[ix2]
  gene_ix <- gene_ix[ix2]
  mu_gene <- mu_gene[ix2]
  FF_gene <- FF_gene[ix2]

  min_vscore <- quantile(Vscores, prob = min_vscore_pctl / 100)

  ix <- ((colSums(E[, gene_ix] >= min_counts) >= min_cells) & (Vscores >= min_vscore))

  if (plot) {

    x_min <- 0.5 * min(mu_gene)
    x_max <- 2 * max(mu_gene)
    xTh <- x_min * exp(log(x_max / x_min) * seq(0, 1, length.out = 100))
    yTh <- (1 + a) * (1 + b) + b * xTh

    # Create a data frame for the line
    line_data <- data.frame(log10_xTh = log10(xTh), log10_yTh = log10(yTh))

    # Create a data frame for points
    points_data <- data.frame(log10_mu_gene = log10(mu_gene)[ix], log10_FF_gene = log10(FF_gene)[ix])

    # Create the ggplot
    g<-ggplot() +
      geom_point(aes(x = log10(mu_gene), y = log10(FF_gene)), col = rgb(0.8, 0.8, 0.8, alpha = 0.3)) +
      geom_point(data = points_data, aes(x = log10_mu_gene, y = log10_FF_gene), col = 'black', alpha = 0.8) +
      geom_line(data = line_data, aes(x = log10_xTh, y = log10_yTh), color = "blue") +
      labs(title = sample_name, x = 'log10(mean)', y = 'log10(Fano factor)')+theme_bw()
    print(g)

  }

  return(gene_ix[ix])
}

#' Calculate variance across the specified axis for a sparse matrix
#'
#' This function computes the variance across the specified axis for a sparse matrix.
#'
#' @param E A sparse matrix where rows represent observations and columns represent features.
#' @param axis An integer specifying the axis along which the variance is calculated (1 for rows, 2 for columns).
#'
#' @return A numeric vector containing the variance values.
#' @importFrom Matrix rowMeans
sparse_var <- function(E, axis = 1) {
  if(axis == 1){
    mean_gene <- rowMeans(E, sparse = TRUE)
  } else {
    mean_gene <- colMeans(E, sparse = TRUE)
  }
  tmp <- E
  tmp@x <- tmp@x^2
  if(axis == 1){
    return(rowMeans(tmp, sparse = TRUE) - mean_gene^2)
  } else {
    return(colMeans(tmp, sparse = TRUE) - mean_gene^2)
  }
}
# sparse_var <- function(E, axis = 1) {
#   mean_gene <- apply(E, axis, mean)
#   tmp <- E
#   tmp@x <- tmp@x^2
#   squared_mean <- apply(tmp, axis, mean)
#   return(squared_mean - mean_gene^2)
# }

#' Multiply each row of a sparse matrix by a scalar
#'
#' This function multiplies each row of a sparse matrix by a scalar value.
#'
#' @param E A sparse matrix where rows represent observations and columns represent features.
#' @param a A scalar value to multiply each row.
#'
#' @return A sparse matrix with each row multiplied by the scalar value.
#' @import Matrix
sparse_multiply <- function(E, a) {
  nrow <- nrow(E)
  w <- sparseMatrix(i = 1:nrow, j = 1:nrow, x = a)
  return(w %*% E)
}

#' Z-score normalize each column of a sparse matrix
#'
#' This function z-score normalizes each column of a sparse matrix.
#'
#' @param E A sparse matrix where rows represent observations and columns represent features.
#' @param gene_mean A vector of mean values for each feature (default is calculated from the matrix).
#' @param gene_stdev A vector of standard deviation values for each feature (default is calculated from the matrix).
#'
#' @return A z-score normalized sparse matrix.
#' @import Matrix
sparse_zscore <- function(E, gene_mean = NULL, gene_stdev = NULL) {
  if (is.null(gene_mean)) {
    gene_mean <- colMeans(E, sparse = TRUE)
  }
  if (is.null(gene_stdev)) {
    gene_stdev <- sqrt(sparse_var(E, axis = 2))
  }
  sm<-matrix(rep(colMeans(E), each = dim(E)[1]), nrow=dim(E)[1])
  return(t(sparse_multiply(t(E - sm), 1 / gene_stdev)))
}

#' Subsample counts in a sparse matrix
#'
#' This function subsamples counts in a sparse matrix based on a given rate.
#'
#' @param E A sparse matrix where rows represent observations and columns represent features.
#' @param rate Subsampling rate for counts (values between 0 and 1).
#' @param original_totals A vector of original total counts for each observation.
#' @param random_seed An integer specifying the random seed for reproducibility.
#'
#' @return A list containing the subsampled sparse matrix (E) and the final downsampling totals.
#' @importFrom Matrix rowSums
#' @importFrom stats rbinom
subsample_counts <- function(E, rate, original_totals, random_seed = 0) {
  if (rate < 1) {
    set.seed(random_seed)
    E@x <- rbinom(length(E@x), round(E@x), rate)
    current_totals <- rowSums(E, sparse = TRUE)
    unsampled_orig_totals <- original_totals - current_totals
    unsampled_downsamp_totals <- rbinom(length(unsampled_orig_totals),
                                        round(unsampled_orig_totals), rate)
    final_downsamp_totals <- current_totals + unsampled_downsamp_totals
  } else {
    final_downsamp_totals <- original_totals
  }
  return(list(E = E, final_downsamp_totals = final_downsamp_totals))
}



#' Return threshold value based on minimum method.
#'
#' The histogram of the input \code{image} is computed if not provided and
#' smoothed until there are only two maxima. Then the minimum in between is
#' the threshold value.
#'
#' Either image or hist must be provided. In case hist is given, the actual
#' histogram of the image is ignored.
#'
#' @param image (M, N[, ...]) numeric vector, optional
#'     Grayscale input image.
#' @param nbins int, optional
#'     Number of bins used to calculate histogram. This value is ignored for
#'     integer arrays.
#' @param max_num_iter int, optional
#'     Maximum number of iterations to smooth the histogram.
#' @param hist array, or 2-tuple of arrays, optional
#'     Histogram to determine the threshold from and a corresponding array
#'     of bin center intensities. Alternatively, only the histogram can be
#'     passed.
#'
#' @return threshold float
#'     Upper threshold value. All pixels with an intensity higher than
#'     this value are assumed to be foreground.
#' @importFrom stats filter
#' @references
#'     C. A. Glasbey, "An analysis of histogram-based thresholding
#'     algorithms," CVGIP: Graphical Models and Image Processing,
#'     vol. 55, pp. 532-537, 1993.
#'     Prewitt, JMS & Mendelsohn, ML (1966), "The analysis of cell
#'     images", Annals of the New York Academy of Sciences 128: 1035-1053
#'     :DOI:`10.1111/j.1749-6632.1965.tb11715.x`


threshold_minimum <- function(image = NULL, nbins = 256, max_num_iter = 10000, hist = NULL) {
  # Helper function to find local maxima indices
  find_local_maxima_idx <- function(hist) {
    maximum_idxs <- c()
    direction <- 1

    for (i in 1:(length(hist) - 1)) {
      if (direction > 0) {
        if (hist[i + 1] < hist[i]) {
          direction <- -1
          maximum_idxs <- c(maximum_idxs, i)
        }
      } else {
        if (hist[i + 1] > hist[i]) {
          direction <- 1
        }
      }
    }

    return(maximum_idxs)
  }

  # Validate image histogram
  result <- validate_image_histogram(image, hist, nbins)
  counts <- result$counts
  bin_centers <- result$bin_centers

  # Initialize smooth_hist
  smooth_hist <- as.numeric(counts)

  # Iterate to smooth the histogram
  for (counter in 1:max_num_iter) {
    smooth_hist <- filter(smooth_hist, rep(1/3, 3), sides = 2)
    maximum_idxs <- find_local_maxima_idx(smooth_hist)

    if (length(maximum_idxs) < 3) {
      break
    }
  }

  # Check if two maxima are found
  if (length(maximum_idxs) != 2) {
    stop('Unable to find two maxima in histogram')
  } else if (counter == max_num_iter) {
    stop('Maximum iteration reached for histogram smoothing')
  }

  # Find the lowest point between the maxima
  threshold_idx <- which.min(smooth_hist[(maximum_idxs[1] + 1):(maximum_idxs[2] + 1)])

  return(bin_centers[maximum_idxs[1] + threshold_idx])
}


#' Helper function to validate image histogram.
#'
#' @param image numeric vector
#'     Grayscale input image.
#' @param hist array, or 2-tuple of arrays
#'     Histogram to determine the threshold from and a corresponding array
#'     of bin center intensities.
#' @param nbins int
#'     Number of bins used to calculate histogram.
#'
#' @return list
#'     A list containing counts and bin_centers.
#'
#' @keywords internal
validate_image_histogram <- function(image, hist, nbins) {
  # Implementation (see previous response)
}

# Helper function to validate image histogram
validate_image_histogram <- function(image, hist, nbins) {
  if (!is.null(hist)) {
    counts <- hist[[1]]
    bin_centers <- hist[[2]]
  } else {
    hist_result <- hist(image, breaks = nbins, plot = FALSE)
    counts <- hist_result$counts
    bin_centers <- hist_result$mids
  }

  return(list(counts = counts, bin_centers = bin_centers))
}

#' @export
print_py<-function(values){
  paste0(paste0(head(values, n =3), collapse=", "), " ... ", paste0(tail(values, n = 3), collapse = " , "))
}
