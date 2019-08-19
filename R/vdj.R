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
#' @importFrom utils read.csv
#'
#' @export
#'
AddClonotype <- function(vdj.dir, scrna) {
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

