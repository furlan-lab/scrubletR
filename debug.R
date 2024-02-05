# debug.R

suppressPackageStartupMessages({
  library(viewmastRust)
  library(Seurat)
  library(SeuratObject)
  library(ggplot2)
  library(scCustomize)
  library(magrittr)
  library(SummarizedExperiment)
})

if(grepl("^gizmo", Sys.info()["nodename"])){
  ROOT_DIR2<-"/fh/fast/furlan_s/grp/data/ddata/BM_data"
} else {
  ROOT_DIR2<-"/Users/sfurlan/Library/CloudStorage/OneDrive-SharedLibraries-FredHutchinsonCancerCenter/Furlan_Lab - General/experiments/patient_marrows/aggr/cds/indy"
}

#query dataset
seuP<-readRDS(file.path(ROOT_DIR2, "220831_WC1.RDS"))
DimPlot_scCustom(seuP)

library(scrubletR)
counts_matrix<-t(seuP@assays$RNA@counts[,1:2000])
counts_matrix<-t(seuP@assays$RNA@counts)



scr<-Scrublet$new(counts_matrix = counts_matrix)
#debug(scr$pipeline_pca)
#debug(scr$scrub_doublets)
#debug(scr$simulate_doublets)
ds<-scr$scrub_doublets()
ds$predicted_doublets

undebug(filter_genes)
seuP$ds<-ds$doublet_scores_obs_

FeaturePlot_scCustom(seuP, features = "ds")


undebug(get_knn_graph)

undebug(sparse_multiply)
debug(scrubletR:::filter_genes)
debug(scr$pipeline_pca)
debug(get_vscores)
debug(scr$scrub_doublets)
debug(scr$simulate_doublets)
undebug(scr$scrub_doublets)
undebug(scr$pipeline_normalize)
undebug(tot_counts_norm)
undebug(sparse_zscore)
debug(scr$pipeline_zscore)
debug(scr$pipeline_normalize_variance)

roxygen2::roxygenise()
devtools::check()

usethis::use_pkgdown_github_pages()
