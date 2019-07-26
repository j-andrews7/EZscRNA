#' Visualize by metadata variables
#'
#' \code{VizDimMetaData} creates Seurat DimPlots for the given meta.data columns
#' using both UMAP, TSNE, and PCA reductions. These plots are automatically 
#' saved in a PDF in the specified output directory.
#'
#' @param scrna Seurat object.
#' @param vars Vector of meta.data column names to plot.
#' @param outdir Path to output directory.
#' @param ... Arguments to be passed to Seurat \code{DimPlot}. \code{cols}, 
#'   \code{reduction}, and \code{group.by} are already defined and will throw
#'   an error if passed.
#'
#' @importFrom Seurat DimPlot
#' @importFrom grDevices dev.off pdf rainbow
#' @import ggplot2
#'
#' @export
#'
VizMetaData <- function(scrna, vars, outdir, ...) {
	message("Plotting meta.data variables.")

	for (i in vars) {
		message("Plotting ", i)
	  # Get all unique elements of variable.
		vars_found <- sort(unique(scrna@meta.data[[as.character(i)]]))
		  
		# Color with rainbow colors.
		cell_colors <- rainbow(length(vars_found), s = 0.6, v = 0.9)

		w <- 10 + (2 * ceiling((length(vars_found) / 12)))

		# Plot.
		pdf(sprintf("%s/UMAP.TSNE.%s.pdf", outdir, i), width = w, height = 5, 
			useDingbats=FALSE)
		p1 <- DimPlot(object = scrna, group.by = as.character(i), 
			cols = cell_colors, reduction = "tsne", ...) + NoLegend()
		p2 <- DimPlot(object = scrna, group.by = as.character(i), 
			cols = cell_colors, reduction = "umap", ...) + 
			theme(legend.text = element_text(size = 7)) +
			labs(color = as.character(i))

		# Finicky garbage to get legends to not overlay plots.
		p2a <- p2 + theme(legend.position = "none")
		legend <- cowplot::get_legend(p2)
		plots <- cowplot::align_plots(p1, p2a, align = 'h', axis = '0')
		leg.w <- ceiling((length(vars_found) / 12)) * 0.4
		c <- cowplot::plot_grid(
		  plots[[1]], plots[[2]], legend,
		  rel_widths = c(1, 1, leg.w),
		  nrow = 1
		)
		print(c)
		dev.off()
	}
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
#' @param ... Arguments to be passed to Seurat \code{DimPlot}. \code{cols}, 
#'   \code{reduction}, and \code{group.by} are already defined and will throw
#'   an error if passed.
#'
#' @importFrom Seurat DimPlot
#' @import ggplot2
#'
#' @export
#'
VizCellType <- function(scrna, outdir, ...) {
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
		reduction = "tsne", ...) + NoLegend()
	p2 <- DimPlot(object = scrna, group.by = "lineage", cols = cell_colors, 
		reduction = "umap", ...) + 
		theme(legend.text = element_text(size = 7)) +
		labs(color = 'Cell Lineage')

	# Finicky garbage to get legends to not overlay plots.
	p2a <- p2 + theme(legend.position = "none")
	legend <- cowplot::get_legend(p2)
	plots <- cowplot::align_plots(p1, p2a, align = 'h', axis = '0')
	leg.w <- ceiling((lineages / 12)) * 0.4
	c <- cowplot::plot_grid(
	  plots[[1]], plots[[2]], legend,
	  rel_widths = c(1, 1, leg.w),
	  nrow = 1
	)
	print(c)
	dev.off()
}


#' Visualize clonotype distributions
#'
#' \code{VizVDJDist} visualizes clonotype distributions for each sample in a
#' seurat object as histograms as well as barcharts comparing clonotype 
#' proportions between them.
#'
#' Cells with no clonotypes are still included in determining clonotype
#' frequencies, but NA is removed from subsequent graphs.
#'
#' @param scrna Seurat object with clonotype data added to metadata with 
#'   \code{AddClonotype}.
#' @param outdir Path to output directory.
#' @param g.by Metadata column to group samples by. If not provided, only
#'   histograms of clonotypes will be saved.
#' @param o.by Vector containing names of members of each group to sort by 
#'   within the group. Ignored if \code{g_by} is NULL. Should contain one 
#'   instance of each potential value in \code{g_by} column if provided.
#' @param n.clono.c Number of top clonotypes to plot for comparison barchart.
#'   Default is 10. Ignored if \code{g.by} is NULL.
#' @param n.clono.g Number of clonotypes to show in group-specific histograms.
#'   All are shown by default.
#'
#' @import ggplot2
#' @import dplyr
#'
#' @export
#'
VizVDJDist <- function(scrna, outdir, g.by = NULL, o.by = NULL, n.clono.c = 10,
	n.clono.g = NULL) {

	message("Visualizing clonotype distributions.")

	if (!is.null(g.by)) {

		# Get clonotype frequencies.
		df <- scrna@meta.data %>% dplyr::count((!!as.symbol(g.by)), cdr3s_aa) %>% 
			dplyr::group_by(!!as.symbol(g.by)) %>% dplyr::mutate(prop = prop.table(n))

		# Sort by frequency.
		df <- df[order(-df$prop),]

		# Remove rows with no TCR data.
		df <- df[!is.na(df$cdr3s_aa),]

		# Get the top n TCR sequences by frequency. Sample agnostic.
		c.df <- df[!duplicated(df$cdr3s_aa),]
		c.df <- c.df[1:as.numeric(n.clono.c),]

		# Get all rows for the the top TCR AA sequences in original df.
		cdr3 <- c.df$cdr3s_aa
		x.df <- df[(df$cdr3s_aa %in% cdr3),]

		# Order 
		index <- 1:length(o.by)
		z = 1
		x.df$ord <- 1
		for (i in o.by) {
		    x.df$ord[x.df[g.by] == i] <- z
		    z = z + 1
		}
		
		pdf(sprintf("%s/Clonotype.Distributions.pdf", outdir), height = 9, 
			width = 9)

		# Plot grouped barchart for top clonotypes.
		p <- ggplot(x.df, aes(x = reorder(cdr3s_aa, -prop), prop, 
			group = ord,  fill = !!as.symbol(g.by))) + 
			geom_col(position = position_dodge(preserve = "single"), 
			colour = "black") + scale_y_continuous(labels = scales::percent) +
			xlab("Clonotype") + ylab("Frequency") + theme_classic() +
			theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), 
			axis.line.x = element_line(size = 1, colour = "black"))
		print(p)

		# Plot individual group frequencies.
		groups <- unique(df[[g.by]])
		for (i in groups) {
			g.spec <- df[which(df[g.by] == i), ]

			# Subset df here for top clonotypes to show in specific groups.
			if (!is.null(n.clono.g)) {
				g.spec <- g.spec[order(-g.spec$prop),]
				g.spec <- g.spec[1:as.numeric(n.clono.g),]
			}

			p <- ggplot(g.spec, aes(x = reorder(cdr3s_aa, -prop), prop)) + geom_col(
				colour = "black") + scale_y_continuous(labels = scales::percent) +
				ggtitle(i) + xlab("Clonotype") + ylab("Frequency") + theme_classic() +
				theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), 
				axis.line.x = element_line(size = 1, colour = "black"))
				print(p)
		}
	} else {
		# Get clonotype frequencies.
		df <- scrna@meta.data %>% dplyr::count((!!as.symbol(g.by)), cdr3s_aa) %>% 
			dplyr::mutate(prop = prop.table(n))

		# Subset df here for top clonotypes to show in specific groups.
		if (!is.null(n.clono.g)) {
			g.spec <- df[order(-df$prop),]
			g.spec <- g.spec[1:as.numeric(n.clono.g),]
		}

		p <- ggplot(g.spec, aes(x = reorder(cdr3s_aa, -prop), prop)) + geom_col(
			colour = "black") + scale_y_continuous(labels = scales::percent) +
			xlab("Clonotype") + ylab("Frequency") + theme_classic() +
			theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), 
			axis.line.x = element_line(size = 1, colour = "black")) 
			print(p)
	}

	dev.off()
}

#' Visualize an annotated marker list
#'
#' \code{VizAnnotatedMarkers} creates Seurat FeaturePlots for lists of 
#' marker genes associated with a given annotation.
#'
#' Plots for each annotation set will be output to a PDF that will be
#' dynamically sized and named based on the name and number of markers for the
#' set. New directories will be created for each set in the output directory.
#'
#' @param scrna Seurat object.
#' @param marker.df Dataframe with the two columns called "Set" and "Marker". 
#'   The "Set" column should contain a cell or process-type (e.g. Tcell, Bcell, 
#'   Exhaustion markers, etc.) while the "Marker" column contains the 
#'   comma-delimited gene symbols associated with it. 
#' @param outdir Path to output directory.
#' @param pt.size Use to adjust point size for plotting FeaturePlots.
#' @param vln Boolean indicating whether to create Seurat VlnPlots for each set.
#'   Splits by cell idents. NULL by default.
#' @param ridge Boolean indicating whether to create Seurat RidgePlots for each
#'   set. Splits by cell idents. NULL by default.
#' @param dot Boolean indicating whether to create Seurat DotPlots for each
#'   set. Splits by cell idents. NULL by default.
#' @param heatmap Boolean indicating whether to create a Seurat Heatmap for each
#'   set. Splits by cell idents. NULL by default.
#'
#' @import Seurat
#' @import ggplot2
#' 
#' @export
#'
VizAnnotatedMarkers <- function(scrna, marker.df, outdir, pt.size = NULL,
	vln = NULL, ridge = NULL, dot = NULL, heatmap = NULL) {

  # Plot individual genes in various classes.
  for (i in unique(marker.df$Set)) {
    # Remove problematic characters from cell classes.
    j <- gsub(" ", "_", i)
    j <- gsub("/", "_", j)
    message("Plotting ", j)
    genes <- trimws(unlist(strsplit(marker.df$Marker[which(marker.df$Set == i)], 
    	",")))
    # Remove genes not found in Seurat object.
    genes <- genes[which(genes %in% rownames(scrna))]
    ng <- length(genes)

    # Move to next set if no genes for current set are in Seurat object.
    if (ng == 0) {
    	message("Skipping ", j, " as no markers are present in Seurat object.")
    	next
    }

    out <- sprintf("%s/%s", outdir, j)
    dir.create(file.path(out), showWarnings = FALSE)
    out.tsne <- sprintf("%s/TSNE.%s.pdf", out, j)
    out.umap <- sprintf("%s/UMAP.%s.pdf", out, j)

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
    fp <- FeaturePlot(object = scrna, features = genes, 
    	cols = c("gray","red"), ncol = 2)
    print(fp)
    dev.off()

    pdf(out.umap, useDingbats = FALSE, height = h, width = w)
    fp <- FeaturePlot(object = scrna, features = genes, 
    	cols = c("gray","red"), ncol = 2, reduction = "umap", pt.size = pt.size)
    print(fp)
    dev.off()

    # Additional plots.
    if (vln) {
    	message("Plotting violin plots.")
    	out.vln <- sprintf("%s/VlnPlots.%s.pdf", out, j)
    	VizVlnPlot(scrna, out.vln, genes)
    }

    if (ridge) {
    	message("Plotting ridge plots.")
    	out.rid <- sprintf("%s/RidgePlots.%s.pdf", out, j)
    	VizRidgePlot(scrna, out.rid, genes)
    }

    if (dot) {
    	message("Plotting dot plots.")
    	out.dot <- sprintf("%s/DotPlots.%s.pdf", out, j)
    	VizDotPlot(scrna, out.dot, genes)
    }

    if (heatmap) {
    	message("Plotting heatmaps.")
    	out.heat <- sprintf("%s/Heatmaps.%s.pdf", out, j)
    	VizHeatmap(scrna, out.heat, genes)
    }

  }
}


#' Visualize scored annotated marker sets
#'
#' \code{VizScoredSets} creates Seurat FeaturePlots for sets of 
#' marker genes associated with a given annotation.
#'
#' Plots for each annotation set will be output to a PDF that will be
#' dynamically named based on the set name. 
#'
#' @param scrna Seurat object.
#' @param sets Vector of 
#' @param outdir Path to output directory.
#' @param pt.size Use to adjust point size for plotting FeaturePlots.
#' @param vln Boolean indicating whether to create Seurat VlnPlots for each set.
#'   Splits by cell idents. NULL by default.
#' @param ridge Boolean indicating whether to create Seurat RidgePlots for each
#'   set. Splits by cell idents. NULL by default.
#' @param dot Boolean indicating whether to create Seurat DotPlots for each
#'   set. Splits by cell idents. NULL by default.
#' @param heatmap Boolean indicating whether to create a Seurat Heatmap for each
#'   set. Splits by cell idents. NULL by default.
#'
#' @import Seurat
#' @import ggplot2
#' 
#' @export
#'
VizScoredSets <- function(scrna, sets, outdir, pt.size = NULL,
	vln = NULL, ridge = NULL, dot = NULL, heatmap = NULL) {

  # Plot individual genes in various classes.
  for (i in sets) {
    # Remove problematic characters from cell classes.
    j <- gsub(" ", "_", i)
    j <- gsub("/", "_", j)
    message("Plotting ", j)
    genes <- trimws(unlist(strsplit(marker.df$Marker[which(marker.df$Set == i)], 
    	",")))
    # Remove genes not found in Seurat object.
    genes <- genes[which(genes %in% rownames(scrna))]
    ng <- length(genes)

    # Move to next set if no genes for current set are in Seurat object.
    if (ng == 0) {
    	message("Skipping ", j, " as no markers are present in Seurat object.")
    	next
    }

    out <- sprintf("%s/%s", outdir, j)
    dir.create(file.path(out), showWarnings = FALSE)
    out.tsne <- sprintf("%s/TSNE.%s.pdf", out, j)
    out.umap <- sprintf("%s/UMAP.%s.pdf", out, j)

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
    fp <- FeaturePlot(object = scrna, features = genes, 
    	cols = c("gray","red"), ncol = 2)
    print(fp)
    dev.off()

    pdf(out.umap, useDingbats = FALSE, height = h, width = w)
    fp <- FeaturePlot(object = scrna, features = genes, 
    	cols = c("gray","red"), ncol = 2, reduction = "umap", pt.size = pt.size)
    print(fp)
    dev.off()

    # Additional plots.
    if (vln) {
    	message("Plotting violin plots.")
    	out.vln <- sprintf("%s/VlnPlots.%s.pdf", out, j)
    	VizVlnPlot(scrna, out.vln, genes)
    }

    if (ridge) {
    	message("Plotting ridge plots.")
    	out.rid <- sprintf("%s/RidgePlots.%s.pdf", out, j)
    	VizRidgePlot(scrna, out.rid, genes)
    }

    if (dot) {
    	message("Plotting dot plots.")
    	out.dot <- sprintf("%s/DotPlots.%s.pdf", out, j)
    	VizDotPlot(scrna, out.dot, genes)
    }

    if (heatmap) {
    	message("Plotting heatmaps.")
    	out.heat <- sprintf("%s/Heatmaps.%s.pdf", out, j)
    	VizHeatmap(scrna, out.heat, genes)
    }

  }
}

# Internal======================================================================

#' Creates VlnPlots for gene set
#'
#' Dynamically generates Seurat VlnPlots for all genes in a given gene set as a
#' PDF.
#'
#' @param scrna Seurat object.
#' @param outfile Path for output PDF.
#' @param genes Vector of genes to plot.
#'
#' @importFrom Seurat VlnPlot
#' @import ggplot2
#'
VizVlnPlot <- function(scrna, outfile, genes) {
	ng <- length(genes)
	if (ng == 1) {
	  w <- 0.3 * length(sort(unique(Idents(scrna))))
	  h <- 4 * ceiling(length(genes) / 3)
	} else if (ng == 2) {
	  w <- 2 * (0.3 * length(sort(unique(Idents(scrna)))))
	  h <- 4 * ceiling(length(genes) / 3)
	} else {
	  w <- 3 * (0.3 * length(sort(unique(Idents(scrna)))))
	  h <- 4 * ceiling(length(genes) / 3)
	}
	pdf(outfile, useDingbats = FALSE, height = h, width = w)
	p <- VlnPlot(scrna, features = genes) + NoLegend()
	print(p)
	dev.off()
}


#' Creates RidgePlots for gene set
#'
#' Dynamically generates Seurat RidgePlots for all genes in a given gene set as
#' a PDF.
#'
#' @param scrna Seurat object.
#' @param outfile Path for output PDF.
#' @param genes Vector of genes to plot.
#'
#' @importFrom Seurat RidgePlot
#' @import ggplot2
#'
VizRidgePlot <- function(scrna, outfile, genes) {
	ng <- length(genes)
	if (ng == 1) {
	    w <- 5
	    h <- 0.2 * length(sort(unique(Idents(scrna))))
	} else if (ng == 2) {
	    w <- 10
	    h <- 0.2 * length(sort(unique(Idents(scrna))))
	} else if (ng == 3) {
	    w <- 15
	    h <- 0.2 * length(sort(unique(Idents(scrna))))
	} else {
	    w <- 15
	    h <- 4 * (0.2 * length(sort(unique(Idents(scrna)))))
	}
	pdf(outfile, useDingbats = FALSE, height = h, width = w)
	p <- RidgePlot(scrna, features = genes) + NoLegend()
	print(p)
	dev.off()
}


#' Creates DotPlot for gene set
#'
#' Dynamically generates a Seurat DotPlot for all genes in a given gene set as
#' a PDF.
#'
#' @param scrna Seurat object.
#' @param outfile Path for output PDF.
#' @param genes Vector of genes to plot.
#'
#' @importFrom Seurat DotPlot
#' @import ggplot2
#'
VizDotPlot <- function(scrna, outfile, genes) {
	ng <- length(genes)
	w <- 5 + (0.3 * ng)
	h <- 0.2 * length(sort(unique(Idents(scrna))))

	pdf(outfile, useDingbats = FALSE, height = h, width = w)
	p <- DotPlot(scrna, features = genes) + theme(axis.text.x = 
		element_text(angle = 45, vjust = 1, hjust = 1))
	print(p)
	dev.off()
}


#' Creates Heatmap for gene set
#'
#' Dynamically generates Seurat Heatmap for all genes in a given gene set as
#' a PDF.
#'
#' @param scrna Seurat object.
#' @param outfile Path for output PDF.
#' @param genes Vector of genes to plot.
#'
#' @importFrom Seurat DotPlot
#' @import ggplot2
#'
VizHeatmap <- function(scrna, outfile, genes) {
	ng <- length(genes)

	h <- 4 + (0.3 * length(genes))
	w <- 3 + (0.4 * length(sort(unique(Idents(scrna)))))

	pdf(outfile, useDingbats = FALSE, height = h, width = w)
	p <- DoHeatmap(subset(scrna, downsample = 100), features = genes, size = 3) + 
		ggtitle("Downsampled to 100 cells")
	print(p)
	p <- DoHeatmap(subset(scrna, downsample = 1000), features = genes, size = 3) + 
		ggtitle("Downsampled to 1000 cells")
	print(p)
	p <- DoHeatmap(scrna, features = genes, size = 3) + 
		ggtitle("No Downsampling")
	print(p)
	dev.off()
}