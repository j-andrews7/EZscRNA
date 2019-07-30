#' Visualize GO/pathway enrichments
#'
#' \code{VizEnrichments} creates barcharts for enrichments returned by
#' \code{\link{RunEnrichr}}. Two plots will be created for each library - one 
#' ranked by adj. p-value, the other by Enrichr's combined score metrics. Plots 
#' will be saved in a PDF named by library.
#'
#' @param enrichments Named list containing enrichment results as returned by
#'   \code{\link{RunEnrichr}}
#' @param outdir Path to the output directory.
#' @param n.terms Number of terms to place on plot. 10 by default.
#' @param remove.insig Boolean indicating whether terms that don't meet the 
#'   \code{adj.p.thresh} should be removed from the plots. TRUE by default.
#' @param adj.p.thresh Value indicating adjusted p-value threshold to filter 
#'   terms. 0.05 by default.
#' @param colors Vector of two colors to use for coloring bars. First will be
#'   low, second will be high. "grey" and "darkred" by default.
#'
#' @import ggplot2
#' @importFrom grDevices dev.off pdf
#' @importFrom scales wrap_format
#'
#' @export
#'
#' @examples
#' genes <- c("CD8A", "GZMA", "GZMK", "GZMB", "CD4", "CD3E", "GNLY")
#' libs <- c("Reactome_2016", "KEGG_2019_HUMAN")
#' terms <- RunEnrichr(genes, libraries = libs)
#'
#' \dontrun{
#' VizEnrichments(enrichments = terms)
#' }
#'
#' @seealso \code{\link{VizEnrichments}} for visualization.
#'
VizEnrichments <- function(enrichments, outdir = "./", 
	n.terms = 10, remove.insig = TRUE, adj.p.thresh = 0.05, 
	colors = c("grey", "darkred")) {

	if (length(colors) > 2) {
		warning("Only two colors can be provided, get that fancy stuff outta here.")
	}

	ind = 1
	for (i in enrichments) {
		lib <- names(enrichments[ind])
		ind <- ind + 1
		# Significance filter.
		if (isTRUE(remove.insig)) {
			i <- i[which(i$Adjusted.P.value <= adj.p.thresh), ]
			if (nrow(i) == 0) {
				message(paste0("Skipping ", lib, " due to no terms meeting the ",
					"significance threshold."))
				next
			}
		}

		message(paste0("Plotting ", lib))

		i$log.Adj.p <- -log10(as.numeric(i$Adjusted.P.value))
		i.p <- i[order(-i$log.Adj.p),]
		i.s <- i[order(-i$Combined.Score),]

		# Limit number of terms.
		if (!is.null(n.terms)) {
			if (nrow(i) > n.terms) {
				i.p <- i.p[1:n.terms, ]
				i.s <- i.s[1:n.terms, ]
			}
		}
		
		p1 <- ggplot(i.p, aes(x = reorder(Term, log.Adj.p),
			log.Adj.p, fill = log.Adj.p)) + geom_col() + coord_flip() +
			cowplot::theme_cowplot(12) + theme(axis.text.y = element_text(size = 7), 
				legend.title = element_text(size=10), 
				legend.text = element_text(size = 10)) + 
			scale_x_discrete(labels = wrap_format(40)) + 
			scale_fill_gradient(low = colors[1], high = colors[2]) + xlab("Term") +
			ylab("-log10(Adjusted p-value)") + ylim(0, NA)

		p2 <- ggplot(i.s, aes(x = reorder(Term, Combined.Score), 
			Combined.Score, fill = log.Adj.p)) + geom_col() + coord_flip() +
			cowplot::theme_cowplot(12) + theme(axis.text.y = element_text(size = 7), 
				legend.title = element_text(size=10), 
				legend.text = element_text(size = 10)) +
			scale_x_discrete(labels = wrap_format(40)) +
			scale_fill_gradient(low = colors[1], high = colors[2]) + xlab("Term") +
			ylab("Combined Score (Enrichr)")

		h <- 1 + 0.4 * nrow(i.p)
		pdf(sprintf("%s/%s.Enrichments.pdf", outdir, lib), height = h, width = 12)
		c <- cowplot::plot_grid(
		  p1, p2, rel_widths = c(1, 1),
		  nrow = 1
		)
		print(c)
		dev.off()
	}
}


#' Visualize by metadata variables
#'
#' \code{VizDimMetaData} creates Seurat DimPlots for the given meta.data columns
#' using UMAP, TSNE, and PCA reductions. These plots are automatically saved in 
#' a PDF in the specified output directory.
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

	# Check for necessary reductions.
	if (is.null(scrna@reductions$tsne)) {
		warning("TSNE not available in reductions, please run RunTSNE() on object.")
		return()
	} else if(is.null(scrna@reductions$umap)) {
		warning("UMAP not available in reductions, please run RunUMAP() on object.")
		return()
	} else if(is.null(scrna@reductions$pca)) {
		warning("PCA not available in reductions, please run RunPCA() on object.")
		return()
	}

	dim.params <- list(...)
	for (i in vars) {
		message("Plotting ", i)
	  # Get all unique elements of variable.
		vars_found <- sort(unique(scrna@meta.data[[as.character(i)]]))
		  
		# Color with rainbow colors.
		cell_colors <- rainbow(length(vars_found), s = 0.6, v = 0.9)

		w <- 15 + (2 * ceiling((length(vars_found) / 12)))

		# Plot.
		pdf(sprintf("%s/PCA.TSNE.UMAP.%s.pdf", outdir, i), width = w, height = 5, 
			useDingbats=FALSE)
		p1 <- do.call(DimPlot, c(scrna, list(group.by = as.character(i), 
			cols = cell_colors, reduction = "pca"), dim.params)) + NoLegend()
		p2 <- do.call(DimPlot, c(scrna, list(group.by = as.character(i), 
			cols = cell_colors, reduction = "tsne"), dim.params)) + NoLegend()
		p3 <- do.call(DimPlot, c(scrna, list(group.by = as.character(i), 
			cols = cell_colors, reduction = "umap"), dim.params)) + 
			theme(legend.text = element_text(size = 7)) +
			labs(color = as.character(i))

		# Finicky garbage to get legends to not overlay plots.
		p3a <- p3 + theme(legend.position = "none")
		legend <- cowplot::get_legend(p3)
		plots <- cowplot::align_plots(p1, p2, p3a, align = 'h', axis = '0')
		leg.w <- ceiling((length(vars_found) / 12)) * 0.4
		c <- cowplot::plot_grid(
		  plots[[1]], plots[[2]], plots[[3]],  legend,
		  rel_widths = c(1, 1, 1, leg.w),
		  nrow = 1
		)
		print(c)
		dev.off()
	}
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
#' \code{VizAnnotatedMarkers} creates Seurat FeaturePlots (and other 
#' visualizations) for lists of marker genes associated with a given annotation.
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
#' @param vln Boolean indicating whether to create Seurat VlnPlots for each set.
#'   Splits by cell idents. NULL by default.
#' @param ridge Boolean indicating whether to create Seurat RidgePlots for each
#'   set. Splits by cell idents. NULL by default.
#' @param dot Boolean indicating whether to create Seurat DotPlots for each
#'   set. Splits by cell idents. NULL by default.
#' @param heatmap Boolean indicating whether to create a Seurat Heatmap for each
#'   set. Splits by cell idents. NULL by default.
#' @param vln.params List of keyword arguments to be passed to Seurat
#'   \code{VlnPlot}. \code{features} are already defined and will throw
#'   an error if passed.
#' @param ridge.params Boolean indicating whether to create Seurat RidgePlots 
#'   for each set. Splits by cell idents. NULL by default.
#' @param dot.params Boolean indicating whether to create Seurat DotPlots for 
#'   each set. Splits by cell idents. NULL by default.
#' @param heatmap.params Boolean indicating whether to create a Seurat Heatmap 
#'   for each set. Splits by cell idents. NULL by default.
#' @param ... Arguments to be passed to Seurat \code{FeaturePlot}. \code{cols}, 
#'   \code{features}, \code{reduction}, and \code{ncol} are already defined and
#'   will throw an error if passed.
#'
#' @import Seurat
#' @import ggplot2
#' 
#' @export
#'
VizAnnotatedMarkers <- function(scrna, marker.df, outdir, vln = NULL, 
	ridge = NULL, dot = NULL, heatmap = NULL, vln.params = NULL, 
	ridge.params = NULL, dot.params = NULL, heatmap.params = NULL, ...) {

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
    	cols = c("gray","red"), ncol = 2, reduction = "tsne", ...)
    print(fp)
    dev.off()

    pdf(out.umap, useDingbats = FALSE, height = h, width = w)
    fp <- FeaturePlot(object = scrna, features = genes, 
    	cols = c("gray","red"), ncol = 2, reduction = "umap", ...)
    print(fp)
    dev.off()

    # Additional plots.
    if (vln) {
    	message("Plotting violin plots.")
    	out.vln <- sprintf("%s/VlnPlots.%s.pdf", out, j)
    	VizVlnPlot(scrna, out.vln, genes, vln.params = vln.params)
    }

    if (ridge) {
    	message("Plotting ridge plots.")
    	out.rid <- sprintf("%s/RidgePlots.%s.pdf", out, j)
    	VizRidgePlot(scrna, out.rid, genes, ridge.params = ridge.params)
    }

    if (dot) {
    	message("Plotting dot plots.")
    	out.dot <- sprintf("%s/DotPlots.%s.pdf", out, j)
    	VizDotPlot(scrna, out.dot, genes, dot.params = dot.params)
    }

    if (heatmap) {
    	message("Plotting heatmaps.")
    	out.heat <- sprintf("%s/Heatmaps.%s.pdf", out, j)
    	VizHeatmap(scrna, out.heat, genes, heatmap.params = heatmap.params)
    }

  }
}


#' Visualize scored gene modules
#'
#' \code{VizScoredSets} creates Seurat FeaturePlots for scored gene sets.
#'
#' Plots for each scored gene set will be output to a PDF that will be
#' dynamically sized and named based on the name and number of markers for the
#' set. New directories will be created for each set in the output directory.
#'
#' @param scrna Seurat object.
#' @param marker.df Dataframe with the two columns called "Set" and "Marker". 
#'   The "Set" column should contain a cell or process-type (e.g. Tcell, Bcell, 
#'   Exhaustion markers, etc.) while the "Marker" column contains the 
#'   comma-delimited gene symbols associated with it. 
#' @param outdir Path to output directory.
#' @param vln Boolean indicating whether to create Seurat VlnPlots for each set.
#'   Splits by cell idents. NULL by default.
#' @param ridge Boolean indicating whether to create Seurat RidgePlots for each
#'   set. Splits by cell idents. NULL by default.
#' @param dot Boolean indicating whether to create Seurat DotPlots for each
#'   set. Splits by cell idents. NULL by default.
#' @param vln.params List of keyword arguments to be passed to Seurat
#'   \code{VlnPlot}. \code{features} are already defined and will throw
#'   an error if passed.
#' @param ridge.params Boolean indicating whether to create Seurat RidgePlots 
#'   for each set. Splits by cell idents. NULL by default.
#' @param dot.params Boolean indicating whether to create Seurat DotPlots for 
#'   each set. Splits by cell idents. NULL by default.
#' @param ... Arguments to be passed to Seurat \code{FeaturePlot}. \code{cols}, 
#'   \code{features}, \code{reduction}, and \code{ncol} are already defined and
#'   will throw an error if passed.
#'
#' @import Seurat
#' @import ggplot2
#' 
#' @export
#'
VizScoredSets <- function(scrna, marker.df, outdir, vln = NULL, 
	ridge = NULL, dot = NULL, vln.params = NULL, 
	ridge.params = NULL, dot.params = NULL, ...) {

  # Plot score for each set.
  sets <- c()
  for (i in unique(marker.df$Set)) {
    # Remove problematic characters from cell classes.
    j <- gsub(" ", "_", i)
    j <- gsub("/", "_", j)
    name <- paste0(j,".Score")

    # Move to next set if no score for current set in Seurat object.
    if (is.null(scrna@meta.data[[name]])) {
    	message("Skipping ", j, " as no score is present in Seurat object.")
    	next
    }

    sets <- append(sets, name)
  }

  ns <- length(sets)
  # Dynamic figure sizing.
  h = 0
  w = 0
  if (ns > 1) {
    w = 11 # two columns
    h = 5 * (floor(ns / 2) + ns %% 2) # number of rows
  } else {
    w = 5 
    h = 5 
  }

  out <- sprintf("%s/ScoredModules", outdir)
  dir.create(file.path(out), showWarnings = FALSE)
  out.tsne <- sprintf("%s/TSNE.ScoredModules.pdf", out)
  out.umap <- sprintf("%s/UMAP.ScoredModules.pdf", out)

  pdf(out.tsne, useDingbats = FALSE, height = h, width = w)
  fp <- FeaturePlot(object = scrna, features = sets, 
  	cols = c("gray","red"), ncol = 2, reduction = "tsne", ...)
  print(fp)
  dev.off()

  pdf(out.umap, useDingbats = FALSE, height = h, width = w)
  fp <- FeaturePlot(object = scrna, features = sets, 
  	cols = c("gray","red"), ncol = 2, reduction = "umap", ...)
  print(fp)
  dev.off()

  # Additional plots.
  if (vln) {
  	message("Plotting violin plots.")
  	out.vln <- sprintf("%s/VlnPlots.ScoredModule.pdf", out)
  	VizVlnPlot(scrna, out.vln, sets, vln.params = vln.params)
  }

  if (ridge) {
  	message("Plotting ridge plots.")
  	out.rid <- sprintf("%s/RidgePlots.ScoredModule.pdf", out)
  	VizRidgePlot(scrna, out.rid, sets, ridge.params = ridge.params)
  }

  if (dot) {
  	message("Plotting dot plots.")
  	out.dot <- sprintf("%s/DotPlots.ScoredModule.pdf", out)
  	VizDotPlot(scrna, out.dot, sets, dot.params = dot.params)
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
#' @param vln.params List of keyword arguments to be passed to Seurat
#'   \code{VlnPlot}. \code{features} are already defined and will throw
#'   an error if passed. NULL by default.
#'
#' @importFrom Seurat VlnPlot
#' @import ggplot2
#'
VizVlnPlot <- function(scrna, outfile, genes, vln.params = NULL) {
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

	# Check for additional kwargs.
	if (!is.null(vln.params)) {
		p <- do.call(VlnPlot, c(scrna, list(features = genes), vln.params)) +
			NoLegend()
	} else {
		p <- VlnPlot(scrna, features = genes) + NoLegend()
	}
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
#' @param ridge.params List of keyword arguments to be passed to Seurat
#'   \code{RidgePlot}. \code{features} are already defined and will throw
#'   an error if passed. NULL by default.
#'
#' @importFrom Seurat RidgePlot
#' @import ggplot2
#'
VizRidgePlot <- function(scrna, outfile, genes, ridge.params = NULL) {
	ng <- length(genes)
	if (ng == 1) {
	    w <- 5
	    h <- 0.3 * length(sort(unique(Idents(scrna))))
	} else if (ng == 2) {
	    w <- 10
	    h <- 0.3 * length(sort(unique(Idents(scrna))))
	} else if (ng == 3) {
	    w <- 15
	    h <- 0.3 * length(sort(unique(Idents(scrna))))
	} else {
	    w <- 15
	    h <- 4 * (0.3 * length(sort(unique(Idents(scrna)))))
	}
	pdf(outfile, useDingbats = FALSE, height = h, width = w)
	# Check for additional kwargs.
	if (!is.null(ridge.params)) {
		p <- do.call(RidgePlot, c(scrna, list(features = genes), ridge.params)) +
			NoLegend()
	} else {
		p <- RidgePlot(scrna, features = genes) + NoLegend()
	}
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
#' @param dot.params List of keyword arguments to be passed to Seurat
#'   \code{DotPlot}. \code{features} are already defined and will throw
#'   an error if passed. NULL by default.
#'
#' @importFrom Seurat DotPlot
#' @import ggplot2
#'
VizDotPlot <- function(scrna, outfile, genes, dot.params = NULL) {
	ng <- length(genes)
	w <- 5 + (0.3 * ng)
	h <- 0.4 * length(sort(unique(Idents(scrna))))

	pdf(outfile, useDingbats = FALSE, height = h, width = w)
	# Check for additional kwargs.
	if (!is.null(dot.params)) {
		p <- do.call(DotPlot, c(scrna, list(features = genes), dot.params)) + 
			theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
	} else {
		p <- DotPlot(scrna, features = genes) + theme(axis.text.x = 
			element_text(angle = 45, vjust = 1, hjust = 1))
	}
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
#' @param heatmap.params List of keyword arguments to be passed to Seurat
#'   \code{DoHeatmap}. \code{features} are already defined and will throw
#'   an error if passed. NULL by default.
#'
#' @importFrom Seurat DotPlot
#' @import ggplot2
#'
VizHeatmap <- function(scrna, outfile, genes, heatmap.params = NULL) {
	ng <- length(genes)

	h <- 4 + (0.3 * length(genes))
	w <- 3 + (0.4 * length(sort(unique(Idents(scrna)))))

	pdf(outfile, useDingbats = FALSE, height = h, width = w)

	# Check for additional kwargs.
	if (!is.null(heatmap.params)) {
		p <- do.call(DoHeatmap, c(subset(scrna, downsample = 100), 
			list(features = genes), heatmap.params)) + 
			ggtitle("Downsampled to 100 cells")
		print(p)
		p <- do.call(DoHeatmap, c(subset(scrna, downsample = 1000), 
			list(features = genes), heatmap.params)) + 
			ggtitle("Downsampled to 1000 cells")
		print(p)
		p <- do.call(DoHeatmap, c(scrna, list (features = genes), heatmap.params)) + 
			ggtitle("Downsampled to 100 cells")
		print(p)
	} else {
		p <- DoHeatmap(subset(scrna, downsample = 100), features = genes, 
			size = 3) + ggtitle("Downsampled to 100 cells")
		print(p)
		p <- DoHeatmap(subset(scrna, downsample = 1000), features = genes, 
			size = 3) + ggtitle("Downsampled to 1000 cells")
		print(p)
		p <- DoHeatmap(scrna, features = genes, size = 3) + 
			ggtitle("No Downsampling")
		print(p)
	}
	dev.off()
}
