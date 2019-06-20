#' Infers and assigns cell type for each cell
#'
#' \code{AssignCellType} performs correlation-based cell inference using a 
#' user-provided reference dataset. It returns either a seurat object with
#' lineage, cell.type, and corr columns added to the metadata, or a dataframe 
#' containing the top three predicted cell types for each cell along with their 
#' correlation values. This dataframe will be saved in the output directory 
#' regardless, along with distribution statistics for inferred cell types.
#'
#' The reference dataset can be from any source, it should just be normalized so
#' that columns (cell types) are comparable. The meat of this code was written
#' by Allegra Petti - the version here has just been made more generic.
#'
#' @param scrna Seurat object.
#' @param dataset Path to tab-delimited table of gene counts. First column must
#'   be gene identifiers of same type as Seurat object. Each subsequent
#'   column should contain counts for a cell type with the column header
#'   denoting the cell type. Replicates should contain an underscore followed 
#'   by the replicate number (e.g. BCell_1, BCell_2, etc).
#' @param outdir Path to output directory.
#' @param assign Boolean indicating whether inferred cell types should actually
#'   be assigned to seurat object or just returned as a table. TRUE by default.
#' @param n_cores Number of cores to use for correlation. Linearly decreases 
#'   computation time.
#' @return If assign is TRUE, returns a seurat object with inferred cell type
#'   information in the metadata. If FALSE, returns a dataframe 
#'
#' @importFrom Seurat AddMetaData
#'
#' @export
#'
AssignCellType <- function(scrna, dataset, outdir, assign = TRUE, n_cores = 1) {
	predictions <- InferCellType(scrna, dataset, outdir, n_cores)

	if (isTRUE(assign)) {

	  print("Incorporating inferred cell types.")

		names(predictions)  <- c("CellBarcode", "Lineage.1", "Lineage.2", 
			"Lineage.3", "Celltype.1", "Celltype.2", "Celltype.3", "Max.1", "Max.2", 
			"Max.3", "Linmed.1", "Linmed.2", "Linmed.3")

		# Most lenient version - take top prediction only.
		final_predictions <- predictions[, c("CellBarcode", "Lineage.1", 
			"Celltype.1", "Max.1")] 

		# Create metadata table and add to seurat object.
		names(final_predictions) <- c("cell", "lineage", "celltype", "corr")
		rownames(final_predictions) <- final_predictions[, 1]
		final_predictions[, 1] <- NULL
		scrna <- AddMetaData(object = scrna, metadata = final_predictions)

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
#'   by the replicate number (e.g. BCell_1, BCell_2, etc).
#' @param outdir Path to output directory.
#' @param n_cores Number of cores to use for correlation. Linearly decreases 
#'   computation time.
#'
#' @return A table containing the top three predicted cell types for each cell.
#'
#' @importFrom foreach foreach
#' @import doParallel
#'
InferCellType <- function(scrna, dataset, outdir, n_cores){
	print(paste("Preparing reference dataset and seurat object for cell type", 
		" inference.", sep = ""))
	x <- PrepRef(scrna, dataset)
	data_sorted <- x$ref
	scrna_subset <- x$sc.sub

	# Set up the dataframe to hold the results.
	cell_barcodes <- colnames(scrna_subset)
	celltype_list <- colnames(data_sorted)

	# Save the lineage types without potential replicate values
	lineages <- gsub("_\\d+", "", names(data_sorted))

	# Set cores for doParallel.
	registerDoParallel(cores = n_cores)

	# Perform correlation for each cell.
	print("Performing cell type inference.")
	results <- foreach(i = 1:length(cell_barcodes), .combine = rbind) %dopar% {
	  if (i%%100 == 0) { 
	    print(paste("Cells: ", i, sep=""))
	  }

	  cell_ind <- which(colnames(scrna_subset) == cell_barcodes[i])
	  profile <- as.vector(scrna_subset[, cell_ind])
	  # List of correlations, one for each column.
	  corrs <- vector(length = ncol(data_sorted), mode = "numeric")
	  k <- 0

	  # Actually do correlation for each column/lineage.
	  for (j in 1:ncol(data_sorted)) {
	  	k <- k+1
	    	corrs[k] <- cor(profile, data_sorted[,j], method = "spearman", 
	    		use = "pairwise.complete.obs")
	  }

	  # Provide the indexes of the corrs vector in descreasing order.
	  order_ind <- order(corrs, decreasing = TRUE)
	  # Get the indices and values for the top three predictions.
	  i1 <- order_ind[1]
	  i2 <- order_ind[2]
	  i3 <- order_ind[3]
	  max1 <- corrs[i1]
	  max2 <- corrs[i2]
	  max3 <- corrs[i3]
	  lineage1 <- as.character(lineages[i1])
	  lineage2 <- as.character(lineages[i2])
	  lineage3 <- as.character(lineages[i3])
	  celltype1 <- as.character(celltype_list[i1])
	  celltype2 <- as.character(celltype_list[i2])
	  celltype3 <- as.character(celltype_list[i3])

	  # Calculate median correlation for top 3 lineages for cell.
	  medians <- vector(length = length(unique(lineages)), mode = "numeric")
	  for (l in 1:length(unique(lineages))) {
	  	lineage_ind <- which(lineages == unique(lineages)[l])
	    medians[l] <- median(corrs[lineage_ind])
	  }
	  sorted_median_ind <- order(medians, decreasing = TRUE)

	  j1 <- sorted_median_ind[1]
	  j2 <- sorted_median_ind[2]
	  j3 <- sorted_median_ind[3]
	  linmed1 <- unique(lineages)[j1]
	  linmed2 <- unique(lineages)[j2]
	  linmed3 <- unique(lineages)[j3]

	  # Stick results in a dataframe.
	  results_temp <- data.frame(cell_barcodes[i], lineage1, lineage2, 
	  	lineage3, celltype1, celltype2, celltype3, max1, max2, max3, 
	    linmed1, linmed2, linmed3)
	}

	# Save output.
	output_table <- table(results$lineage1)
	output_table_pct <- 100 * (table(results$lineage1) / 
		sum(table(results$lineage1)))
	write.table(output_table, file = sprintf("%s/TopLineage.Distribution.txt", 
		outdir), quote = FALSE, sep = "\t", row.names = FALSE)
	write.table(output_table_pct, 
		file = sprintf("%s/TopLineage.Distribution.Pct.txt", 
		outdir), quote = FALSE, sep = "\t", row.names = FALSE)
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
#'
#' @return A list of two dataframes: "ref" - sorted genes from the reference 
#'   dataset also found in the seurat object, and "sc.sub" - sorted genes from 
#'   the seurat object also found in the reference dataset.
#'
#' @import dplyr
#' @import Seurat
#'
PrepRef <- function(scrna, dataset){
	# Get both the genes and cells for the normalized seurat object.
	genes <- rownames(x = scrna)
	cells <- Cells(scrna)

	# Load reference dataset.
	data_all <- read.table(dataset, sep = "\t", row.names = 1, header = TRUE,
		fill = TRUE, quote = '', check.names = FALSE)
	# Only keep genes found in the seurat object as well.
	data_unsorted <- data_all[which(rownames(data_all) %in% genes),]

	# Extract signature genes from normalized seurat object.
	g_ind = which(genes %in% rownames(data_unsorted))
	scrna_subset = FetchData(object=scrna, vars=genes[g_ind])
	# Transpose so rows are genes, columns are cells. 
	scrna_subset = t(scrna_subset)

	# Sort the data so the genes are in the same order as in the subsetted 
	# seurat object.
	new_indices = match(rownames(scrna_subset),rownames(data_unsorted))
	data_sorted <- data_unsorted[new_indices,]

	return(list("ref" = data_sorted, "sc.sub" = scrna_subset))
}
