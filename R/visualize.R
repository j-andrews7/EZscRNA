#' Visualize an annotated marker list
#'
#' \code{VizAnnotatedMarkers} creates Seurat \code{FeaturePlot} for lists of 
#' marker genes associated with a given annotation.
#'
#' Plots for each annotation set will be output to a PDF that will be
#' dynamically sized and named based on the name and number of markers for the
#' set.
#'
#' @param scrna Seurat object.
#' @param marker.file Tab-delimited file with the first two columns named "Set"
#'   and "Marker". The "Set" column should contain a cell or process-type (e.g. 
#'   Tcell, Bcell, Exhaustion markers, etc.) while the "Marker" column contains
#'   the gene symbol associated with it. Many rows can be added to visualize 
#'   many "Markers" for each "Set". 
#' @param outdir Path to output directory.
#' 
#' @export
#'
VizAnnotatedMarkers <- function(scrna, marker.file, outdir) {
	# Read and filter the gene lists.
	message("Parsing and filtering gene lists.")
	glists.raw <- read.table(marker.file, sep = "\t", row.names = NULL, 
		header = TRUE, as.is = TRUE)
  glists <- glists.raw[which(glists.raw$Marker %in% rownames(scrna)), ] 

  # Plot individual genes in various classes.
  message("Plotting gene sets.")
  for (i in unique(glists$Set)) {
    gplotlist <- list()
    # Remove problematic characters from cell classes.
    j <- gsub(" ", "_", i)
    j <- gsub("/", "_", j)
    message("Plotting ", j)
    genesToPlot <- glists$Marker[which(glists$Set == i)]
    ng <- length(genesToPlot) 
    out.tsne <- sprintf("%s/TSNE.%s.pdf", outdir, j)
    out.umap <- sprintf("%s/UMAP.%s.pdf", outdir, j)

    # Dynamic figure sizing.
    h = 0
    w = 0
    if (ng > 1) {
      w = 10 # two columns
      h = 5 * (floor(ng / 2) + ng %% 2) # number of rows
    } else {
      w = 5 
      h = 5 
    }
    pdf(out.tsne, useDingbats = FALSE, height = h, width = w)
    fp <- FeaturePlot(object = scrna, features = genesToPlot, 
    	cols = c("gray","red"), ncol = 2, reduction = "tsne")
    print(fp)
    dev.off()

    pdf(out.umap, useDingbats = FALSE, height = h, width = w)
    fp <- FeaturePlot(object = scrna, features = genesToPlot, 
    	cols = c("gray","red"), ncol = 2, reduction = "umap")
    print(fp)
    dev.off()
  }
}

#' Visualize by metadata variables
#'
#'
#'
#'
#'
#'
#'
#'
#' @export
#'
VizMetaData <- function(scrna, vars, outdir) {

}

#' Visualize inferred cell types
#'
#' \code{VizCellType} creates Seurat DimPlots for the cell types and lineages
#' using both UMAP and TSNE reductions. These plots are automatically saved in
#' a PDF in the specified output directory.
#'
#' @param scrna Seurat object with base.lineage, lineage, and celltype columns
#'   added to the meta.data (as done by \code{AssignCellType}).
#' @param outdir Path to output directory.
#'
#'
#' @importFrom Seurat DimPlot
#' @import cowplot
#'
#' @export
#'
VizCellType <- function(scrna, outdir) {
  message("Plotting cell types.")
  # All found lineages.
	lineage_found <- sort(unique(scrna@meta.data$lineage))
	  
	# Color with rainbow colors.
	cell_colors <- rainbow(length(lineage_found), s = 0.6, v = 0.9)

	w <- 10 + (2 * ceiling((length(lineage_found) / 12)))

	# Plot.
	pdf(sprintf("%s/UMAP.TSNE.celltype.pdf", outdir), width = w, height = 5, 
		useDingbats=FALSE)
	p1 <- DimPlot(object = scrna, group.by = "lineage", cols = cell_colors, 
		pt.size = 0.3, reduction = "tsne") + NoLegend()
	p2 <- DimPlot(object = scrna, group.by = "lineage", cols = cell_colors, 
		pt.size = 0.3, reduction = "umap") + 
		theme(legend.text = element_text(size = 7)) +
		labs(color='Cell Lineage')

	# Finicky garbage to get legends to not overlay plots.
	p2a <- p2 + theme(legend.position = "none")
	legend <- get_legend(p2)
	plots <- align_plots(p1, p2a, align = 'h', axis = '0')
	leg.w <- ceiling((lineages / 12)) * 0.4
	c <- plot_grid(
	  plots[[1]], plots[[2]], legend,
	  labels = c("A", "B"),
	  rel_widths = c(1, 1, leg.w),
	  nrow = 1
	)
	print(c)
	dev.off()
}