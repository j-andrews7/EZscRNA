#' Score annotated marker lists
#'
#' \code{ScoreAnnotatedMarkers} scores \linkS4class{Seurat} objects for lists of 
#' marker genes associated with a given annotation. Scores are stored in
#' \code{meta.data} columns named after each set.
#'
#' @param scrna \linkS4class{Seurat} object.
#' @param marker.df Dataframe with two columns named "Set" and "Marker". 
#'   The "Set" column should contain a cell or process-type (e.g. Tcell, Bcell, 
#'   Exhaustion markers, etc.) while the "Marker" column contains the 
#'   comma-delimited gene symbols associated with it.
#' @return \linkS4class{Seurat} object with scores for each "Set" added to 
#'   \code{meta.data}.
#'
#' @importFrom Seurat AddModuleScore
#' 
#' @export
#'
#' @seealso \code{\link[Seurat]{AddModuleScore}} for scoring info and 
#'   \code{\link{VizAnnotatedMarkers}} for convenient visualization options.
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
    	), assay = "RNA")

    # Remove random number from column header.
    # No idea why, but this does not work with the "[[]]" accessor.
    colnames(scrna@meta.data)[colnames(scrna@meta.data) == paste0(j, ".Score1"
    	)] <- paste0(j, ".Score")
  }

  return(scrna)
}