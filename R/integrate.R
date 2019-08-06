#' Integrate Seurat objects simply and easily
#'
#' \code{SimpleIntegration} integrates Seurat objects using default Seurat 
#' parameters that should work for most cases. It can take either a list of 
#' Seurat objects or a meta.data variable to split a single Seurat object by
#' prior to integration. 
#'
#' Performs \code{SCTransform} on each Seurat object individually, per Seurat's
#' vignette for best practices
#'
#' @param scrnas Either a list of Seurat objects or a single Seurat object, the
#'   latter of which requires \code{split.by} to be provided as well.
#' @param split.by String containing name of meta.data column to be used to 
#'   split Seurat object into multiple samples. 
#' @param skip.SCT Boolean indicating whether the \code{SCTransform} call for
#'   each sample should be skipped or not. This should only be set to TRUE if
#'   the samples if the samples have already undergone SCTransform. 
#' @param vars.to.regress Vector of meta.data variables to be regressed out
#'   during \code{SCTransform}. 
#' @param n.features Number of anchor features to use. 3000 by default.
#' @return An integrated Seurat object with technical variation/batch effects
#'   removed from the individual samples.
#'
#' @import Seurat
#' @import future
#'
#' @export
#'
SimpleIntegration <- function(scrnas, split.by = NULL, skip.SCT = FALSE, 
	vars.to.regress = NULL, n.features = 3000) {

	# Parameter checking.
	if (length(scrnas) > 1 & !isTRUE(skip.SCT)) {
		for (i in 1:length(scrnas)) {
			scrnas[[i]] <- SCTransform(scrnas[[i]], verbose = FALSE, 
				vars.to.regress = vars.to.regress)
		}
	} else if (is.null(split.by)) {
		stop(paste0("split.by must be provided if a list of Seurat objects is,",
			"not provided as input."))
	} else {
		scrnas <- SplitObject(scrnas, split.by = split.by)
		for (i in 1:length(scrnas)) {
			scrnas[[i]] <- SCTransform(scrnas[[i]], verbose = FALSE, 
				vars.to.regress = vars.to.regress)
		}
	}

	# Actual integration.
	anch.features <- SelectIntegrationFeatures(object.list = scrnas, 
		nfeatures = n.features)
	scrnas <- PrepSCTIntegration(object.list = scrnas, 
		anchor.features = anch.features, verbose = FALSE)

	anchors <- FindIntegrationAnchors(object.list = scrnas, 
		normalization.method = "SCT", anchor.features = anch.features, 
		verbose = FALSE)
	scrna <- IntegrateData(anchorset = anchors, 
		normalization.method = "SCT", verbose = FALSE)

	return(scrna)
}