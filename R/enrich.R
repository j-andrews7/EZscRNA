#' Perform enrichment analyses for a gene list
#'
#' \code{RunEnrichr} submits gene lists to the enrichr web server for gene 
#' ontology/pathway enrichment analyses and collects the results as a named list
#' for each library used.
#'
#' Potential libraries can be viewed with \code{listEnrichrDbs} from the
#' \code{enrichR} package.
#'
#' The enrichr web server can be found at 
#' \url{https://amp.pharm.mssm.edu/Enrichr/}. If you use this function, you 
#' should be sure to cite the original authors.
#'
#' @param genes A vector of gene symbols that will be tested.
#' @param libraries A vector of libraries to test the genes against.
#' @param outdir Path to output directory. If specified, tables for each library
#'   will be saved.
#' @return A named list of enrichment results for each library.
#'
#' @importFrom enrichR enrichr
#' @importFrom utils write.table
#'
#' @export
#' 
#' @references \href{https://www.ncbi.nlm.nih.gov/pubmed/27141961}{Kuleshov MV, 
#'   Jones MR, Rouillard AD, et al. Enrichr: a comprehensive gene set enrichment 
#'   analysis web server 2016 update. Nucleic Acids Research. 2016; gkw377}
#'
#' @examples
#' genes <- c("CD8A", "CD4", "FOXP3", "CTLA4", "GNLY")
#' libs <- c("Reactome_2016", "KEGG_2019_HUMAN")
#' enrichments <- RunEnrichr(genes, libraries = libs)
#'
#' @seealso \code{\link{VizEnrichments}} for visualization.
#'
RunEnrichr <- function(genes, libraries = c("GO_Molecular_Function_2018", 
	"GO_Cellular_Component_2018", "GO_Biological_Process_2018", "KEGG_2019_Human",
	"Reactome_2016", "BioCarta_2016", "Panther_2016"), outdir = NULL) {

	# Required or package will throw an error saying the website is not responding
	# due to poor checking.
	options(base_address = "http://amp.pharm.mssm.edu/Enrichr/")
	options(enrichRLive = TRUE)

	# Run enrichments.
	message("Submitting sets to the Enrichr server. Please be patient, this can",
		" sometimes take a few minutes if the server is under heavy load or down.")
	results <- enrichr(genes, libraries)
	res.names <- names(results)

	if (!is.null(outdir)) {
		y = 1
		for (i in results) {
			df <- as.data.frame(i, sep = "\t")
			df$Old.P.value <- NULL
			df$Old.Adjusted.P.value <- NULL
			write.table(df, file = sprintf("%s/%s.Results.txt", outdir, res.names[y]),
				sep = "\t", quote = FALSE, row.names = FALSE)
			y <- y + 1
		}
	}

	return(results)
}