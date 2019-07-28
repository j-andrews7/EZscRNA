#' Normalize, scale, and regress out wanted variation
#'
#' \code{ClusterDEG} runs \code{SCTransform} on a Seurat object, followed by 
#' PCA, TSNE, UMAP, and clustering. Also finds marker genes and saves the output 
#' as a table along with a heatmap of the top 10 upregulated genes in each 
#' cluster.
#'
#' If multiple \code{res} values are given, a table and heatmap will be make for
#' each, along with saving the cluster numbers for each in their own meta.data
#' columns.
#'
#' @param scrna Seurat object.
#' @param outdir Path to output directory.
#' @param npcs Number of principle components to use for UMAP and clustering.
#'   Default is 30, as \code{SCTransform} tends to do better with more.
#' @param res Numeric value denoting resolution to use for clustering. Default
#'   is 0.8 (Seurat default). Increasing this value will speed up clustering,
#'   but may decrease the numbers of distinct clusters. Values of 0.5-3 are 
#'   sensible. Multiple values may be entered as a vector - resulting clusters
#'   will be added as a meta.data column named "Cluster_res_npcs" where "res"
#'   and "npcs" will be the values for those arguments, respectively. The
#'   clusters derived from the last value in the list will be set as the
#'   default Ident for cells. 
#' @param skip.sct Boolean indicating whether to skip \code{SCTransform}. FALSE
#'   by default. Set to TRUE if \code{SimpleIntegration} was used to integrate
#'   the Seurat object.
#' @param min.dist Number that controls how tighly the embedding is allowed to
#'   compress points together in \code{RunUMAP}. Increasing may be beneficial
#'   for large datasets. Default is 0.3 (Seurat default).
#' @param n.neighbors Integer that determines the number of neighboring points
#'   used in local approximations of manifold structure in \code{RunUMAP}.
#'   Altering it may be beneficial for large datasets (though it isn't stated
#'   how it should be changed). Values of 5-50 are considered sensical.
#'   Default is 30 (Seurat default).
#' @param regress Vector of metadata variables to regress during data scaling.
#'   Must match column headers in meta.data. 
#' @param ccpca Boolean to indicate whether PCA using only cell cycle genes
#'   should be done. If so, it will saved as a reduction named "cc". FALSE by
#'   default.
#' @param test Denotes which DE test to use for marker finding. Options are: 
#'   "wilcox" (default), "bimod", "roc", "t", "negbinom", "poisson", "LR", 
#'   "MAST", "DESeq2".
#' @param logfc.thresh Value that limits DE testing to genes that show, 
#'   on average, at least X-fold difference (log-scale) between two groups of 
#'   cells. Increasing speeds up function at cost of potentially missing weaker 
#'   differences. Default is 0.25 (Seurat default).
#' @param min.pct Value that limits DE testing to genes detected in a minimum
#'   fraction of cells in either population. Default is 0.1 (Seurat default).
#' @return A Seurat object with normalized, scaled counts and assigned clusters.
#'   If \code{ccpca = TRUE}, an additional PCA reduction named "cc" will also be 
#'   present.
#'
#' @import Seurat
#' @import sctransform
#' @importFrom grDevices dev.off pdf
#' @import ggplot2
#' @importFrom utils write.table
#'
#' @export
#'
ClusterDEG <- function(scrna, outdir, npcs = 30, res = 0.8, skip.sct = FALSE, 
									min.dist = 0.3, n.neighbors = 30, regress = NULL, 
									ccpca = FALSE, test = "wilcox", 
									logfc.thresh = 0.25, min.pct = 0.1) {

  # Run sctransform & regress out any specified variables.
  if(!skip.sct) {
  	message("Scaling & normalizing data with SCTransform.")
  	scrna <- SCTransform(scrna, vars.to.regress = regress, 
  		return.only.var.genes = FALSE )
	}

  # Run PCA using just cell cycle genes if indicated. Saved as "cc".
  if (isTRUE(ccpca)) {
		message("Performing PCA on cell cycle genes.")
		s.genes <- cc.genes$s.genes
		g2m.genes <- cc.genes$g2m.genes
		scrna <- RunPCA(scrna, npcs = npcs, features = c(s.genes, g2m.genes),
		  reduction.name = "cc")
  }

  message("Performing PCA/UMAP/TSNE on variable features.")
  scrna <- RunPCA(scrna, npcs = npcs)
  scrna <- RunUMAP(scrna, dims = 1:npcs)
  scrna <- RunTSNE(scrna, dims = 1:npcs)

  message("Performing clustering on variable features.")
  scrna <- FindNeighbors(scrna, dims = 1:npcs)

  for (i in res) {
  	message(paste0("Finding cluster markers using ", res, "resolution."))
  	scrna <- FindClusters(scrna, resolution = res)
  	scrna[[sprintf("Clusters.%.1fRes.%dPC", res, npcs)]] <- Idents(object = 
  		scrna) 
  
	  markers <- FindAllMarkers(scrna, assay = "RNA", 
			logfc.threshold = logfc.thresh, min.pct = min.pct, test.use = test)
		
	    # Save markers as table.
	  message("Saving cluster markers.")
	  write.table(markers, file = sprintf("%s/Cluster.Markers.%.1fRes.%dPC.txt", 
	  	outdir, res, npcs), sep = "\t", quote = FALSE, row.names = FALSE)

	  message("Visualizing cluster markers.")
	  top10up <- markers %>% dplyr::group_by(cluster) %>% 
			dplyr::top_n(n = 10, wt = avg_logFC) %>% dplyr::filter(avg_logFC > 0)
	  # Make marker heatmap. Dynamic sizing.
	  h <- 4 + (0.3 * length(top10up$gene))
		w <- 3 + (0.4 * length(sort(unique(Idents(scrna)))))
	  pdf(sprintf("%s/Top10.UpMarkers.Cluster.%.1fRes.%dPC.Heatmap.pdf", outdir, 
	  	res, npcs), height = h, width = w)
	  DoHeatmap(scrna, features = top10up$gene)
	  dev.off()
	}

  return(scrna)
}
