#' Infers and assigns cell type for each cell
#'
#' \code{AssignCellType} performs correlation-based cell inference using a 
#' user-provided reference dataset. It returns either a seurat object with
#' lineage, cell.type, corr, and (potentially) base.lineage columns added to the 
#' metadata, or a dataframe containing the top three predicted cell types for 
#' each cell along with their correlation values. This dataframe will be saved 
#' in the output directory regardless, along with distribution statistics for 
#' inferred cell types and lineages if set.
#'
#' The reference dataset can be from any source, it should just be normalized so
#' that columns (cell types) are comparable. The majority of this code was 
#' written by Allegra Petti - the version here has just been made more generic.
#'
#' @param scrna Seurat object.
#' @param dataset Path to tab-delimited table of gene counts. First column must
#'   be gene identifiers of same type as Seurat object. Each subsequent
#'   column should contain counts for a cell type with the column header
#'   denoting the cell type. Replicates should contain an underscore followed 
#'   by the replicate number (e.g. BCell_1, BCell_2, etc).
#' @param outdir Path to output directory.
#' @param lineage Boolean indicating whether or not column names are formatted
#'   as lineage followed by more granular specifications (e.g. Tcell.CD4_1,
#'   Tcell.CD8_1, etc). If TRUE, will assign the lineage (everything before the)
#'   first '.' in the column name to a metadata column called "base.lineage".
#'   Columns not containing a '.' will use the entire name, though replicate
#'   indicators will be removed (e.g. NK_1 will just be NK).
#' @param assign Boolean indicating whether inferred cell types should actually
#'   be assigned to seurat object or just returned as a table. TRUE by default.
#' @param n.cores Number of cores to use for correlation. Linearly decreases 
#'   computation time.
#' @return If assign is TRUE, returns a seurat object with inferred cell type
#'   information in the metadata. If FALSE, returns a dataframe of the 
#'   inferred cell type/lineage information instead.
#'
#' @importFrom Seurat AddMetaData
#'
#' @export
#'
AssignCellType <- function(scrna, dataset, outdir, lineage = FALSE, 
													 assign = TRUE, n.cores = 1) {
	predictions <- InferCellType(scrna, dataset, outdir, lineage, n.cores)

	if (isTRUE(assign)) {

	  message("Incorporating inferred cell types.")

	  if (lineage) {
			names(predictions)  <- c("CellBarcode", "Lineage.1", "Lineage.2", 
				"Lineage.3", "Celltype.1", "Celltype.2", "Celltype.3", "Max.1", "Max.2", 
				"Max.3", "Linmed.1", "Linmed.2", "Linmed.3", "base.lineage")

			# Most lenient version - take top prediction only.
			final.predictions <- predictions[, c("CellBarcode", "Lineage.1", 
				"Celltype.1", "Max.1", "base.lineage")] 

			# Create metadata table and add to seurat object.
			names(final.predictions) <- c("cell", "lineage", "celltype", "corr", 
				"base.lineage")
		} else {
			# Same thing, but with no base lineage.
			names(predictions)  <- c("CellBarcode", "Lineage.1", "Lineage.2", 
				"Lineage.3", "Celltype.1", "Celltype.2", "Celltype.3", "Max.1", "Max.2", 
				"Max.3", "Linmed.1", "Linmed.2", "Linmed.3")

			final.predictions <- predictions[, c("CellBarcode", "Lineage.1", 
				"Celltype.1", "Max.1")] 

			names(final.predictions) <- c("cell", "lineage", "celltype", "corr")
		}

		rownames(final.predictions) <- final.predictions[, 1]
		final.predictions[, 1] <- NULL
		scrna <- AddMetaData(object = scrna, metadata = final.predictions)

	} else {
		scrna <- predictions
	}

	return(scrna)
}


# Internal======================================================================

#' Infer cell type using reference dataset
#'
#' \code{InferType} utilizes a reference dataset to perform correlations of each
#' cell in a seurat object with each sample in the reference dataset. It retuns
#' a table containing the most likely cell type for each cell based on 
#' correlation values from the reference dataset. 
#'
#' The reference dataset can be from any source, it should just be normalized so
#' that columns (cell types) are comparable.
#'
#' @param scrna Seurat object.
#' @param dataset Path to tab-delimited table of gene counts. First column must
#'   be gene identifiers that match those of the Seurat object. Each subsequent
#'   column should contain counts for a cell type with the column header
#'   denoting the cell type. Replicates should contain an underscore followed 
#'   by the replicate number (e.g. Bcell_1, Bcell_2, etc). 
#' @param outdir Path to output directory.
#' @param lineage Boolean indicating whether or not column names are formatted
#'   as lineage followed by more granular specifications (e.g. Tcell.CD4_1,
#'   Tcell.CD8_1, etc). If TRUE, will assign the lineage (everything before the)
#'   first '.' in the column name to a metadata column called "base.lineage".
#'   Columns not containing a '.' will use the entire name, though replicate
#'   indicators will be removed (e.g. NK_1 will just be NK).
#' @param n.cores Number of cores to use for correlation. Linearly decreases 
#'   computation time.
#' @return A table containing the top three predicted cell types for each cell.
#'
#' @importFrom foreach foreach
#' @import doParallel
#' @importFrom stats cor median reorder
#' @importFrom utils read.csv read.table write.table
#'
InferCellType <- function(scrna, dataset, outdir, lineage = FALSE, n.cores = 1){
	message(paste("Preparing reference dataset and seurat object for cell type", 
		" inference.", sep = ""))
	x <- PrepRef(scrna, dataset)
	data.sorted <- x$ref
	scrna.subset <- x$sc.sub

	# Set up the dataframe to hold the results.
	cell.barcodes <- colnames(scrna.subset)
	celltype.list <- colnames(data.sorted)

	# Save the lineage types without potential replicate values
	lineages <- gsub("_\\d+", "", names(data.sorted))

	# Set cores for doParallel.
	registerDoParallel(cores = n.cores)

	# Perform correlation for each cell.
	message("Performing cell type inference.")
	results <- foreach(i = 1:length(cell.barcodes), .combine = rbind) %dopar% {
	  if (i%%100 == 0) { 
	    message(paste("Cells: ", i, sep=""))
	  }

	  cell.ind <- which(colnames(scrna.subset) == cell.barcodes[i])
	  profile <- as.vector(scrna.subset[, cell.ind])
	  # List of correlations, one for each column.
	  corrs <- vector(length = ncol(data.sorted), mode = "numeric")
	  k <- 0

	  # Actually do correlation for each column/lineage.
	  for (j in 1:ncol(data.sorted)) {
	  	k <- k+1
	    	corrs[k] <- cor(profile, data.sorted[,j], method = "spearman", 
	    		use = "pairwise.complete.obs")
	  }

	  # Provide the indexes of the corrs vector in descreasing order.
	  order.ind <- order(corrs, decreasing = TRUE)
	  # Get the indices and values for the top three predictions.
	  i1 <- order.ind[1]
	  i2 <- order.ind[2]
	  i3 <- order.ind[3]
	  max1 <- corrs[i1]
	  max2 <- corrs[i2]
	  max3 <- corrs[i3]
	  lineage1 <- as.character(lineages[i1])
	  lineage2 <- as.character(lineages[i2])
	  lineage3 <- as.character(lineages[i3])
	  celltype1 <- as.character(celltype.list[i1])
	  celltype2 <- as.character(celltype.list[i2])
	  celltype3 <- as.character(celltype.list[i3])

	  # Calculate median correlation for top 3 lineages for cell.
	  medians <- vector(length = length(unique(lineages)), mode = "numeric")
	  for (l in 1:length(unique(lineages))) {
	  	lineage.ind <- which(lineages == unique(lineages)[l])
	    medians[l] <- median(corrs[lineage.ind])
	  }
	  sorted.median.ind <- order(medians, decreasing = TRUE)

	  j1 <- sorted.median.ind[1]
	  j2 <- sorted.median.ind[2]
	  j3 <- sorted.median.ind[3]
	  linmed1 <- unique(lineages)[j1]
	  linmed2 <- unique(lineages)[j2]
	  linmed3 <- unique(lineages)[j3]

	  # Stick results in a dataframe.
	  results.temp <- data.frame(cell.barcodes[i], lineage1, lineage2, 
	  	lineage3, celltype1, celltype2, celltype3, max1, max2, max3, 
	    linmed1, linmed2, linmed3)
	}

	# Save output.
	output.table <- table(results$lineage1)
	output.table.pct <- 100 * (table(results$lineage1) / 
		sum(table(results$lineage1)))
	write.table(output.table, file = sprintf("%s/TopCelltype.Distribution.txt", 
		outdir), quote = FALSE, sep = "\t", row.names = FALSE)
	write.table(output.table.pct, 
		file = sprintf("%s/TopCelltype.Distribution.Pct.txt", 
		outdir), quote = FALSE, sep = "\t", row.names = FALSE)

	# Get base lineages.
	if(lineage) {

		results$base.lineage <- lapply(results$lineage1, SplitIt)
		results$base.lineage <- as.character(results$base.lineage)

		# Save output.
		output.table <- table(results$base.lineage)
		output.table.pct <- 100 * (table(results$base.lineage) / 
			sum(table(results$base.lineage)))
		write.table(output.table, file = sprintf("%s/TopLineage.Distribution.txt", 
			outdir), quote = FALSE, sep = "\t", row.names = FALSE)
		write.table(output.table.pct, 
			file = sprintf("%s/TopLineage.Distribution.Pct.txt", 
			outdir), quote = FALSE, sep = "\t", row.names = FALSE)
	}

	write.table(results, file=sprintf("%s/CellType.Predictions.txt", outdir), 
		sep = "\t", quote = FALSE, row.names = FALSE)

	return(results)
}


#' Load and prepare reference dataset for cell type inference
#'
#' \code{PrepRef} processes a reference dataset of normalized gene counts for
#' correlation analysis with a seurat object. Removes unnecessary genes from the
#' reference dataset and matches the row ordering of the seurat object.
#'
#' @param scrna Seurat object.
#' @param dataset Path to tab-delimited table of gene counts. First column must
#'   be gene identifiers that match those of the Seurat object. Each subsequent
#'   column should contain counts for a cell type with the column header
#'   denoting the cell type. Replicates contain an underscore followed by the
#'   replicate number (e.g. BCell_1, BCell_2, etc).
#' @return A list of two dataframes: "ref" - sorted genes from the reference 
#'   dataset also found in the seurat object, and "sc.sub" - sorted genes from 
#'   the seurat object also found in the reference dataset.
#'
#' @importFrom Seurat Cells FetchData
#'
PrepRef <- function(scrna, dataset){
	# Get both the genes and cells for the normalized seurat object.
	genes <- rownames(x = scrna)
	cells <- Cells(scrna)

	# Load reference dataset.
	data.all <- read.table(dataset, sep = "\t", row.names = 1, header = TRUE,
		fill = TRUE, quote = '', check.names = FALSE)
	# Only keep genes found in the seurat object as well.
	data.unsorted <- data.all[which(rownames(data.all) %in% genes), ]

	# Extract signature genes from normalized seurat object.
	g.ind = which(genes %in% rownames(data.unsorted))
	scrna.subset = FetchData(object = scrna, vars = genes[g.ind])
	# Transpose so rows are genes, columns are cells. 
	scrna.subset = t(scrna.subset)

	# Sort the data so the genes are in the same order as in the subsetted 
	# seurat object.
	new.indices = match(rownames(scrna.subset), rownames(data.unsorted))
	data.sorted <- data.unsorted[new.indices,]

	return(list("ref" = data.sorted, "sc.sub" = scrna.subset))
}


#' Split lineages to get base lineages
#'
#' For use with \code{lapply}. Returns only the first element of resulting list.
#'
#' @param x String to be split.
#'
SplitIt <- function(x) {
  unlist(strsplit(as.character(x), ".", 
		fixed = TRUE))[c(TRUE, FALSE, FALSE, FALSE)]
}
