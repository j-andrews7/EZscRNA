#' Normalize, scale, and regress out wanted variation
#'
#' \code{ClusterDEG} runs \code{\link[Seurat]{SCTransform}} on a 
#' \linkS4class{Seurat} object, followed by \code{\link[Seurat]{RunPCA}}, 
#' \code{\link[Seurat]{RunTSNE}}, \code{\link[Seurat]{RunUMAP}}, and 
#' clustering. Also finds marker genes for each cluster and saves the output as
#' a table along with a heatmap of the top 10 upregulated genes in each cluster.
#'
#' @details
#' If multiple \code{res} values are given, a table and heatmap will be made for
#' each, along with saving the clusters for each in their own \code{meta.data} 
#' columns. 
#'
#' Heatmaps created by \code{ClusterDEG} have each identity class downsampled 
#' to a max of 100 cells - this makes smaller clusters much more visible.
#'
#' @param scrna \linkS4class{Seurat} object.
#' @param outdir Path to output directory.
#' @param npcs Number of principle components to use for UMAP and clustering.
#' @param res Numeric value denoting resolution to use for clustering. 
#'   Higher values generally mean fewer clusters. Values of 0.5-3 are sensible. 
#'   Multiple values may be entered as a vector - resulting clusters will be 
#'   added as a \code{meta.data} column named \code{Cluster_res_npcs} where 
#'   \code{res} and \code{npcs} will be the values for those arguments, 
#'   respectively. The clusters derived from the last value in the list will be 
#'   set as the default Ident for cells and stored in \code{meta.data} under
#'   'seurat_clusters' in addition to the aforementioned format.
#' @param mnn Boolean indicating whether \code{scrna} was integrated with
#'   \code{method="MNN"} via \code{\link{SimpleIntegration}}. If so, must be set
#'   to \code{TRUE} or unintegrated PCA embeddings will be used for
#'   dimensionality reduction and clustering.
#' @param skip.sct Boolean indicating whether to skip 
#'   \code{\link[Seurat]{SCTransform}}. Set to TRUE if 
#'   \code{\link{SimpleIntegration}} was used to integrate the 
#'   \code{Seurat} object.
#' @param min.dist Number that controls how tighly the embedding is allowed to
#'   compress points together in \code{\link[Seurat]{RunUMAP}}. Increasing may 
#'   be beneficial for large datasets.
#' @param n.neighbors Integer that determines the number of neighboring points
#'   used in local approximations of manifold structure in 
#'   \code{\link[Seurat]{RunUMAP}}. Values of 5-50 
#'   are considered sensical. Larger values preserve more global structure while 
#'   detailed local structure is lost.
#' @param regress Character vector of \code{meta.data} variables to regress 
#'   during data scaling.
#' @param ccpca Boolean to indicate whether PCA using only cell cycle genes
#'   should be done. If so, it will saved as a reduction named "cc". This is
#'   useful to compare to prior PCAs using the cell cycle genes if cell cycle
#'   scores were regressed out via \code{regress}.
#' @param test String indication which DE test to use for marker finding. 
#'   Options are: 
#'     "wilcox", "bimod", "roc", "t", "negbinom", "poisson", "LR", 
#'     "MAST", "DESeq2". See \code{\link[Seurat]{FindAllMarkers}}.
#' @param logfc.thresh Value that limits DE testing to genes that show, 
#'   on average, at least X-fold difference (log-scale) between two groups of 
#'   cells. Increasing speeds up function at cost of potentially missing weaker 
#'   differences. 
#' @param min.pct Value that limits DE testing to genes detected in a minimum
#'   fraction of cells in either population.
#' @return A \linkS4class{Seurat} object with normalized, scaled counts and 
#'   assigned clusters. If \code{ccpca = TRUE}, an additional PCA reduction 
#'   named "cc" will also be present.
#'
#' @import Seurat
#' @importFrom grDevices dev.off pdf
#' @importFrom utils write.table
#' @import sctransform
#'
#' @export
#'
#' @examples
#' \dontrun{
#' library(Seurat)
#' scrna <- ClusterDEG(pbmc_small)
#' }
#'
#' \dontrun{
#' # Multiple clustering resolutions
#' scrna <- ClusterDEG(pbmc_small, res = c(0.8, 1, 1.2))
#' }
#'
#' @author Jared Andrews
#'
ClusterDEG <- function(scrna, outdir = ".", npcs = 30, res = 0.8, mnn = FALSE,
  skip.sct = FALSE, min.dist = 0.3, n.neighbors = 30, regress = NULL, 
	ccpca = FALSE, test = "wilcox", logfc.thresh = 0.25, min.pct = 0.1) {

  # Run sctransform & regress out any specified variables.
  if(!skip.sct) {
  	message("Scaling & normalizing data with SCTransform.")
  	scrna <- SCTransform(scrna, vars.to.regress = regress)
	}

  # Run PCA using just cell cycle genes if indicated. Saved as "cc".
  if (isTRUE(ccpca)) {
		message("Performing PCA on cell cycle genes.")
		s.genes <- cc.genes$s.genes
		g2m.genes <- cc.genes$g2m.genes
		scrna <- RunPCA(scrna, npcs = npcs, features = c(s.genes, g2m.genes),
		  reduction.name = "cc")
  }

  # Change to "mnn" reduction for fastMNN integrated objects.
  message("Performing PCA/UMAP/TSNE on variable features.")
  if (!mnn) {
    reduc <- "pca"
  } else {
    reduc <- "mnn"
  }
  message("Using ", reduc, " reduction.")

  if (!mnn) {
  scrna <- RunPCA(scrna, npcs = npcs)
  }
  scrna <- RunTSNE(scrna, dims = 1:npcs, reduction = reduc)
  scrna <- RunUMAP(scrna, dims = 1:npcs, n.neighbors = n.neighbors, min.dist =
    min.dist, reduction = reduc, umap.method = "umap-learn", metric = "correlation")

  message("Performing clustering on variable features.")
  scrna <- FindNeighbors(scrna, dims = 1:npcs, reduction = reduc)

  # Clustering/marker finding for multiple resolutions.
  for (i in res) {
  	message(paste0("Finding cluster markers using ", i, " resolution."))
  	scrna <- FindClusters(scrna, resolution = i)
  	scrna[[sprintf("Clusters.%.1fRes.%dPC.%s", i, npcs, reduc)]] <- 
      Idents(scrna) 
  
	  markers <- FindAllMarkers(scrna, assay = "RNA", 
			logfc.threshold = logfc.thresh, min.pct = min.pct, test.use = test)
		
	   # Save markers as table.
	  message("Saving cluster markers.")
	  write.table(markers, 
      file = sprintf("%s/Cluster.Markers.%.1fRes.%dPC.%s.txt", 
	  	outdir, i, npcs, reduc), sep = "\t", quote = FALSE, row.names = FALSE)

	  message("Visualizing cluster markers.")
	  top10up <- markers %>% dplyr::group_by(cluster) %>% 
			dplyr::top_n(n = 10, wt = avg_logFC) %>% dplyr::filter(avg_logFC > 0)
	  # Make marker heatmap. Dynamic sizing.
	  h <- 4 + (0.15 * length(top10up$gene))
		w <- 3 + (0.4 * length(sort(unique(Idents(scrna)))))
	  pdf(sprintf("%s/Top10.UpMarkers.Cluster.%.1fRes.%dPC.%s.Downsample.Heatmap.pdf", 
      outdir, i, npcs, reduc), height = h, width = w)
	  p <- DoHeatmap(subset(scrna, downsample = 100), features = top10up$gene, 
      assay = "RNA")
    print(p)
	  dev.off()
	}

  return(scrna)
}
