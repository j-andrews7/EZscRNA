#' Add 10X clonotype data to seurat object
#'
#' \code{AddClonotype} adds clonotype info from matched 10X VDJ sequencing to 
#' the metadata of a given Seurat object.
#'
#' This function is admittedly rough and will be rewritten in the future. It
#' does not include specific VDJ genes for each cell, rather just using the
#' final amino acid sequence for inter-sample comparison.
#'
#' @param vdj.dir String containing path to TCR (VDJ) directory.
#' @param scrna Seurat object.
#' @return Seurat object with clonotype data (clonotype_id and cdr3s_aa) 
#'   added to the metadata for each cell.
#'
#' @importFrom Seurat AddMetaData
#'
#' @export
#'
AddClonotype <- function(vdj.dir, scrna){
	message("Adding clonotype data.")
	vdj <- read.csv(sprintf("%s/filtered_contig_annotations.csv", vdj.dir))

	# Remove the -1 at the end of each barcode added by cellranger aggr fuction.
	# Subsets so only the first line of each barcode is kept,
	# as each entry for given barcode will have same clonotype id.
	vdj$barcode <- gsub("-1", "", vdj$barcode)
	vdj <- vdj[!duplicated(vdj$barcode), ]

	# Only keep the barcode and clonotype columns. 
	# We'll get additional clonotype info from the clonotype table.
	vdj <- vdj[, c("barcode", "raw_clonotype_id")]
	names(vdj)[names(vdj) == "raw_clonotype_id"] <- "clonotype_id"

	# Clonotype-centric info.
	clono <- read.csv(sprintf("%s/clonotypes.csv", vdj.dir))

	# Slap the AA sequences onto our original table by clonotype_id.
	vdj <- merge(vdj, clono[, c("clonotype_id", "cdr3s_aa")])

	# Reorder so barcodes are first column and set them as rownames.
	vdj <- vdj[, c(2,1,3)]
	rownames(vdj) <- vdj[, 1]
	vdj[, 1] <- NULL

	# Add to the Seurat object's metadata.
	clono.seurat <- AddMetaData(object = scrna, metadata = vdj)
	message("Clonotype data added.")
		
	return(clono.seurat)
}

#' Visualize clonotype distributions
#'
#' \code{VizVDJDist} visualizes clonotype distributions for each sample in a
#' seurat object as histograms as well as barcharts comparing clonotype 
#' proportions between them.
#'
#' Rows with 
#'
#' @param scrna Seurat object with clonotype data added to metadata with 
#'   \code{AddClonotype}.
#' @param outdir Path to output directory.
#' @param g.by Metadata column to group samples by. If not provided, only
#'   histograms of clonotypes will be saved.
#' @param o.by Vector containing names of members of each group to sort by 
#'   within the group. Ignored if g_by is NULL. Should contain one instance of
#'   each potential value in g_by column if provided.
#' @param n.clono.c Number of top clonotypes to plot for comparison barchart.
#'   Default is 10. Ignored if \code(group_by) is NULL.
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
			xlab("Clonotype") + ylab("Frequency") +
			theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
			theme_classic()
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
				ggtitle(i) + xlab("Clonotype") + ylab("Frequency") +
				theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
				theme_classic()
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
			xlab("Clonotype") + ylab("Frequency") +
			theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
			theme_classic()
			print(p)
	}

	dev.off()
}
