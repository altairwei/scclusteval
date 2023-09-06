
#' Calculate Silhouette width from PCA space for each cell after clustering
#' This is calculated from Seurat object
#' @param object A Seurat object with Idents set to cluster ids (factors)
#' @param dims default 1:50  dimension to use in the PCA space to calculate
#' @param reduction Name of reduction to calculate Silhouette on.
#' eucledian distance
#'
#' @return a dataframe with silhouette width for each cell. see also \code{\link[cluster]{silhouette}}
#' @export
#'
#' @examples
#' CalculateSilhouette(pbmc_small, dims = 1:15)
#'
CalculateSilhouette <- function(object, dims = 1:50, reduction = "pca"){
  embed <- Seurat::Embeddings(object, reduction = reduction)
  if (length(dims) > ncol(embed)) {
    stop("please specify PCA dims smaller than calculated")
  }

  cell_distance <- dist(embed[, dims])
  cell_cluster <- Seurat::Idents(object)
  silhouette_score <- cluster::silhouette(as.integer(cell_cluster), cell_distance)
  silhouette_score <- tibble::tibble(cluster = cell_cluster,
                                     width = silhouette_score[, 3],
                                     cell = names(cell_cluster))
  return(silhouette_score)
}
