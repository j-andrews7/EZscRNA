#' Integrate Seurat objects simply and easily
#'
#' \code{SimpleIntegration} integrates \linkS4class{Seurat} objects using 
#' either Seurat integration functions or 
#' \code{RunFastMNN}, a wrapper from \code{SeuratWrappers} around
#' \code{fastMNN} from batchelor. This is useful for removing batch effects. 
#' It can take either a list of \linkS4class{Seurat}  objects or a 
#' \code{meta.data} variable to split a single \linkS4class{Seurat} object by 
#' prior to integration. 
#'
#' @details
#' Performs \code{\link[Seurat]{SCTransform}} on each Seurat object 
#' individually prior to integration. 
#'
#' If \code{method="MNN"} is used, \code{SeuratWrappers} must be installed.
#'
#' @param scrnas Either a list of \linkS4class{Seurat} objects or a single 
#'   \linkS4class{Seurat} object, the latter of which requires \code{split.by} 
#'   to be provided as well.
#' @param split.by String containing name of meta.data column to be used to 
#'   split Seurat object into multiple samples. 
#' @param method String indicating method to be used. "Seurat" will use Seurat's
#'   internal functions for SCT integration. "MNN" will use
#'   \code{RunFastMNN}.
#' @param vars.to.regress Vector of meta.data variables to be regressed out
#'   during \code{\link[Seurat]{SCTransform}}. 
#' @param n.features Number of anchor features to use with 
#'   \code{\link[Seurat]{SelectIntegrationFeatures}}. 
#' @return An integrated \linkS4class{Seurat} object with technical 
#'   variation/batch effects removed from the individual samples.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' library(Seurat)
#' # Will error due to too few cells.
#' pbmc.int <- SimpleIntegration(pbmc_small, split.by = "groups")
#' }
#'
#' @author Jared Andrews
#'
SimpleIntegration <- function(scrnas, split.by = NULL, 
  method = c("Seurat", "MNN"), vars.to.regress = NULL, n.features = 3000) {

  # Arg check.
  method <- match.arg(method)

	# Seurat object splitting checking.
	if (length(scrnas) > 1) {
		for (i in 1:length(scrnas)) {
			scrnas[[i]] <- Seurat::SCTransform(scrnas[[i]], verbose = FALSE, 
				vars.to.regress = vars.to.regress, return.only.var.genes = FALSE)
		}
	} else if (is.null(split.by)) {
		stop(paste0("split.by must be provided if a list of Seurat objects is,",
			"not provided as input."))
	} else {
		scrnas <- Seurat::SplitObject(scrnas, split.by = split.by)
		for (i in 1:length(scrnas)) {
			scrnas[[i]] <- Seurat::SCTransform(scrnas[[i]], verbose = FALSE, 
				vars.to.regress = vars.to.regress, return.only.var.genes = FALSE)
		}
	}

  if (method == "MNN") {
    scrna <- RunFastMNN(object.list = scrnas)
  } else if (method == "Seurat") {
    # Actual integration.
    anch.features <- Seurat::SelectIntegrationFeatures(object.list = scrnas, 
      nfeatures = n.features)
    scrnas <- Seurat::PrepSCTIntegration(object.list = scrnas, 
      anchor.features = anch.features, verbose = FALSE)

    anchors <- Seurat::FindIntegrationAnchors(object.list = scrnas, 
      normalization.method = "SCT", anchor.features = anch.features, 
      verbose = FALSE)
    scrna <- Seurat::IntegrateData(anchorset = anchors, 
      normalization.method = "SCT", verbose = FALSE)
  }

	return(scrna)
}
