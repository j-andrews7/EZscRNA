#' Score annotated marker lists
#'
#' \code{ScoreAnnotatedMarkers} creates Seurat FeaturePlots for lists of 
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
#' @return Seurat object with scores for each "Set" added to meta.data.
#'
#' @importFrom Seurat AddModuleScore
#' 
#' @export
#'
ScoreAnnotatedMarkers <- function(scrna, marker.df) {

  # Plot individual genes in various classes.
  for (i in unique(marker.df$Set)) {

    # Remove problematic characters from cell classes.
    j <- gsub(" ", "_", i)
    j <- gsub("/", "_", j)
    message("Scoring ", j)
    genes <- trimws(unlist(strsplit(marker.df$Marker[which(marker.df$Set == i)], 
    	",")))
    # Remove genes not found in Seurat object.
    genes <- list(genes[which(genes %in% rownames(scrna))])
    ng <- length(genes[[1]])

    # Move to next set if no genes for current set are in Seurat object.
    if (ng == 0) {
    	message("Skipping ", j, " as no markers are present in Seurat object.")
    	next
    }

    # Score set.
    scrna <- AddModuleScore(scrna, features = genes, name = paste0(j, ".Score"
    	))

    # Remove random number from column header, annoys me.
    colnames(scrna[[]])[colnames(scrna[[]]) == paste0(j, ".Score1"
    	)] <- paste0(j, ".Score")
  }

  return(scrna)
}