#' Visualize an annotated marker list
#'
#'
#'
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
VizAnnotatedMarkers <- function(scrna, marker_file, outdir) {
	# gene lists
	glists.raw <- read.table(gene.lists, sep = ",", row.names = NULL, 
		header = TRUE, as.is = TRUE)
  glists <- glists.raw[which(glists.raw$Name %in% rownames(scrna)), ] # filtered genelists
                                         # plot individual genes in various classes:
  for (i in unique(glists$List)) {
    gplotlist=list()
    j=gsub(" ","_",i)
    j=gsub("/","_",j)
    genesToPlot = glists$Name[which(glists$List == i)]
    ng = length(genesToPlot) # number of genes
    outfile = sprintf("%s/tsneplot.%s.%s.%s.pdf", output.genes, j, control, 
    	date)
    print(outfile)
    h = 0
    w = 0
    if (ng > 1) {
      w = 10 # two columns
      h = 5 * (floor(ng / 2) + ng %% 2) # number of rows
    } else {
      w = 5 # two columns
      h = 5 # number of rows
      }
    pdf(outfile, useDingbats = FALSE, height = h, width = w)
    fp <- FeaturePlot(object = scrna, features = genesToPlot, 
    	cols = c("gray","red"), ncol=2, reduction = "tsne")
    print(fp)
    dev.off()
  }
}

#' Visualize a list of genes
#'
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
VizGeneList <- function(scrna, genelist, outdir) {

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
#'
#'
#'
#'
#'
#'
#'
#'
#'
#' @importFrom Seurat DimPlot
#'
#' @export
#'
VizCellType <- function(scrna, outdir) {
	# # Make t-SNE plot colored according to lineage prediction.
	n <- length(lineage_list)

	# All found lineages.
	lineage_found <- sort(unique(scrna@meta.data$lineage))
	  
	# Color with rainbow colors.
	rainbow_colors <- rainbow(n, s = 0.6, v = 0.9)
	cell_colors <- rainbow_colors[which(lineage_list %in% lineage_found)] 

	# Plot.
	pdf(sprintf("%s/UMAP.TSNE.celltype.pdf", outdir), width = 12, height = 5, 
		useDingbats=FALSE)
	p1 <- DimPlot(object = scrna, group.by = "lineage", cols = cell_colors, 
		pt.size = 0.3, reduction = "tsne")
	p2 <- DimPlot(object = scrna, group.by = "lineage", cols = cell_colors, 
		pt.size = 0.3, reduction = "umap")
	print(CombinePlots(plots = c(p1, p2), legend = "right"))
	dev.off()
}