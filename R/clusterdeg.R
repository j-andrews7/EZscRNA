#' Normalize, scale, and regress out wanted variation
#'
#' \code{SCT} runs \code{SCTransform} on a Seurat object, followed by PCA, UMAP,
#' and clustering. Produces PCA and UMAP dimplots based on user-provided list
#' of metadata columns to use for grouping. Also finds marker genes and
#' saves the output as a table along with a heatmap of the top 10 upregulated
#' genes in each cluster.
#'
#'
#' @param scrna Seurat object.
#' @param outdir Path to output directory.
#' @param npcs Number of principle components to use for UMAP and clustering.
#'   Default is 50, as \code{SCTransform} tends to do better with more.
#' @param res Numeric value denoting resolution to use for clustering. Default
#'   is 0.8 (Seurat default). Increasing this value will speed up clustering,
#'   but may decrease the numbers of distinct clusters. Values of 0.5-3 are 
#'   sensible.
#' @param min_dist Number that controls how tighly the embedding is allowed to
#'   compress points together in \code{RunUMAP}. Increasing may be beneficial
#'   for large datasets. Default is 0.3 (Seurat default).
#' @param n_neighbors Integer that determines the number of neighboring points
#'   used in local approximations of manifold structure in \code{RunUMAP}.
#'   Altering it may be beneficial for large datasets (though it isn't stated
#'   how it should be changed). Values of 5-50 are considered sensical.
#'   Default is 30 (Seurat default).
#' @param regress Vector of metadata variables to regress during data scaling.
#'   Must match column headers in metadata. 
#' @param groups Vector of metadata variables to use for grouping in UMAP and 
#'   PCA dimplots. Must match column headers in metadata. 
#' @param groups_pca Vector of boolean values to determine if \code{DimPlots} 
#'   PCA reductions should also be created for each group. NULL by default (
#'   only UMAP reductions will be shown for each group). If provided, must be
#'   same length as groups parameter.
#' @param groups_label Vector of boolean values to determine if \code{DimPlots} 
#'   should show labels for each group. NULL by default (labels will not be 
#'   shown). If provided, must be same length as groups parameter.
#' @param groups_legend Vector of boolean values to determine if \code{DimPlots} 
#'   should show legends for each group. NULL by default (legends will be 
#'   shown). If provided, must be same length as groups parameter.
#' @param ccpca Boolean to indicate whether PCA using only cell cycle genes
#'   should be performed and plotted by sample identity and Phase.
#' @param use_augment Boolean to indicate whether \code{AugmentPlot} should be
#'   used so points aren't saved as vector graphics while axes, labels, etc are.
#'   Useful if the plot is going to be edited in Illustrator. True by default.
#' @param test Denotes which DE test to use for marker finding. Options are: 
#'   "wilcox" (default), "bimod", "roc", "t", "negbinom", "poisson", "LR", 
#'   "MAST", "DESeq2".
#' @param logfc_thresh Value that limits DE testing to genes that show, 
#'   on average, at least X-fold difference (log-scale) between two groups of 
#'   cells. Increasing speeds up function at cost of potentially missing weaker 
#'   signals. Default is 0.25 (Seurat default).
#' @param min_pct Value that limits DE testing to genes detected in a minimum
#'   fraction of cells in either population. Default is 0.1 (Seurat default).
#'
#' @return A Seurat object with normalized, scaled counts and assigned clusters.
#'
#' @import Seurat
#' @import sctransform
#' @import dplyr
#'
#' @export
#'
RunSCT <- function(scrna, outdir, npcs = 50, res = 0.8, min_dist = 0.3, 
									n_neighbors = 30, regress = NULL, groups = NULL, 
									groups_pca = NULL, groups_label = NULL, groups_legend = NULL,
									ccpca = FALSE, use_augment = TRUE, test = "wilcox", 
									logfc_thresh = 0.25, min_pct = 0.1) {

  # Run sctransform & regress out any specified variables.
  scrna <- SCTransform(scrna, vars.to.regress = regress)

  # Run PCA using just cell cycle genes if indicated. Saved as "cc".
  if (isTRUE(ccpca)) {
		print("Performing PCA on cell cycle genes.")
		s_genes <- cc.genes$s.genes
		g2m_genes <- cc.genes$g2m.genes
		scrna <- RunPCA(scrna, npcs = npcs, features = c(s_genes, g2m_genes),
		  reduction.name = "cc")
		pdf(sprintf("%s/CellCycle.PCA.Regression.pdf", outdir), height = 5, 
		width = 9)
		p <- DimPlot(scrna, group.by = "orig.ident", reduction = "cc", 
			pt.size = 0.3) + ggplot2::theme(legend.position = "right")
		if (isTRUE(use_augment)) {
		  p <- AugmentPlot(plot = p) # Saves points as a png.
		}
		print(p)
		p <- DimPlot(scrna, group.by = "Phase", reduction = "cc", pt.size = 0.3) +
			ggplot2::theme(legend.position = "right")
		if (isTRUE(use_augment)) {
		  p <- AugmentPlot(plot = p) # Saves points as a png.
		}
		print(p)
		dev.off()
  }

  print("Performing PCA/UMAP on variable features.")
  scrna <- RunPCA(scrna, npcs = npcs)
  scrna <- RunUMAP(scrna, dims = 1:npcs)

  print("Performing clustering on variable features.")
  scrna <- FindNeighbors(scrna, dims = 1:npcs)
  scrna <- FindClusters(scrna, resolution = res)
  scrna <- StashIdent(object = scrna, 
	save.name = sprintf("ClusterNames_%.1f_%dPC", res, npcs))

  print("Plotting PCA/UMAP dimplots.")
  pdf(sprintf("%s/EDA.UMAP.PCA.pdf", outdir), height = 5, width = 7)

  # Some preeeetty gross logic here.
  if (!is.null(groups)) {
		for (i in 1:length(groups)) {
		  var = groups[i]
		  label = FALSE
		  p2 = NULL
			if (!is.null(groups_label)) {
			  label = groups_label[i]
			}
			# Yikes.
			if (!is.null(groups_legend)) {

			  if (!isTRUE(groups_legend[i])) {
					p <- DimPlot(scrna, label = label, group.by = var, pt.size = 0.3) +
				  	NoLegend() + ggtitle(label = sprintf("UMAP - %s", var))
					if (!is.null(groups_pca)) {
				  	if (isTRUE(groups_pca[i])) {
						p2 <- DimPlot(scrna, label = label, group.by = var,
					  	reduction = "pca", pt.size = 0.3) + NoLegend() +
					  	ggtitle(label = sprintf("PCA - %s", var))
				  	}
					}
			  } else {
			    p <- DimPlot(scrna, label = label, group.by = var, pt.size = 0.3) +
				  ggtitle(label = sprintf("UMAP - %s", var))
			    if (!is.null(groups_pca)) {
				  	if (isTRUE(groups_pca[i])) {
				    	p2 <- DimPlot(scrna, label = label, group.by = var,
					  	reduction = "pca", pt.size = 0.3) +
					  	ggtitle(label = sprintf("PCA - %s", var))
				  	}
			    }
			  }

			} else {
			  p <- DimPlot(scrna, label = label, group.by = var, pt.size = 0.3) +
				ggtitle(label = sprintf("UMAP - %s", var))
				if (!is.null(groups_pca)) {
			  	if (isTRUE(groups_pca[i])) {
						p2 <- DimPlot(scrna, label = label, group.by = var,
				  		reduction = "pca", pt.size = 0.3) +
				  		ggtitle(label = sprintf("PCA - %s", var))
			  	}
				}
		  }

		  if (isTRUE(use_augment)) {
		    p <- AugmentPlot(plot = p)
		  }
		  print(p)

		  if (!is.null(p2)) {
		    if (isTRUE(use_augment)) {
			  	p2 <- AugmentPlot(plot = p2)
		    }
		    print(p2)
		  }
    }
  }

  p <- DimPlot(scrna, label = TRUE, pt.size = 0.3) + NoLegend() + 
  	ggtitle(label = "UMAP")
	if (isTRUE(use_augment)) {
		p <- AugmentPlot(plot = p)
	}
	print(p)
  

  dev.off()

  print("Finding cluster markers.")
  markers <- FindAllMarkers(scrna, assay = "RNA", 
		logfc.threshold = logfc_thresh, min.pct = min_pct, test.use = test)
	
  # Make marker heatmap.
  pdf(sprintf("%s/Top10.UpMarkers.Cluster.Heatmap.pdf", outdir), height=24, 
	width=38)
  top10up <- markers %>% group_by(cluster) %>% 
		top_n(n = 10, wt = avg_logFC) %>% filter(avg_logFC > 0)
  DoHeatmap(scrna, features = top10up$gene, lines.width=5)
  dev.off()

  # Save markers as table.
  write.table(markers, file=sprintf("%s/Cluster.Markers.txt", outdir), 
		sep="\t", quote=FALSE, row.names=FALSE)

  return(scrna)
}
