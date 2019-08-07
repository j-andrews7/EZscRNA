#' Infers and assigns cell type for each cell using SingleR
#'
#' \code{AssignCellType} performs correlation-based cell inference using a 
#' reference dataset via \code{SingleR}.
#'
#' @details
#' Reference datasets available for use by this function via \code{SingleR} 
#' include:
#' \describe{
#'   \item{"hpca"}{Human Primary Cell Atlas (HPCA): A collection of Gene 
#'     Expression Omnibus (GEO datasets), which contains 713 microarray 
#'     samples classified to 38 main cell types and further annotated to 169
#'     subtypes.}
#'   \item{"blueprint_encode"}{Blueprint + ENCODE datasets: Blueprint 
#'     Epigenomics, 144 RNA-seq pure immune samples annotated to 28 cell 
#'     types. ENCODE, 115 RNA-seq pure stroma and immune samples annotated 
#'     to 17 cell types. Altogether, 259 samples with 43 cell types.}
#'   \item{"immgen"}{Immunological Genome Project (ImmGen): 830 microarray 
#'     samples, which we classified to 20 main cell types and further
#'     annotated to 253 subtypes.}
#'   \item{"mouse.rnaseq"}{A dataset of 358 mouse RNA-seq samples annotated to
#'     28 cell types. This dataset was collected, processed, and shared 
#'     courtesy of Bérénice Benayoun. This data set is especially useful for 
#'     brain-related samples.}
#' }
#'
#' @param scrna Seurat object.
#' @param refset String indicating which reference dataset to download and use
#'   as the training set. 
#' @param labels If \code{TRUE}, use broad lineage-based labels 
#'   \code{main_types} for each reference sample. This is faster and can be more 
#'   informative depending on how specific your data is. If \code{FALSE}, more
#'   specific labels \code{types} will be used.
#' @param method String specifying whether annotation should be applied to 
#'   each single cell \code{method = "single"} or aggregated into cluster-level 
#'   profiles \code{method = "cluster"} prior to annotation.
#' @param clusters String defining the \code{meta.data} column in \code{scrna}
#'   that specify cluster identities for each cell. Only used if \code{method 
#'   = "cluster"}.
#' @param assign Boolean indicating whether inferred cell types should actually
#'   be assigned to \code{Seurat} object or just returned as a \code{dataframe}. 
#' @param n.cores Number of cores to use for correlation. Linearly decreases 
#'   computation time.
#' @param ... Additional arguments to be passed to 
#'   \code{\link[SingleR]{SingleR}}.
#' @return If \code{assign = TRUE}, returns a \code{Seurat} object with inferred
#'   cell type information in the \code{meta.data} slot named in 
#'   \code{refset.labels} format. If \code{method = "cluster"}, the resulting
#'   \code{meta.data} column will be named in \code{clusters.refset.labels} 
#'   format. If \code{assign = FALSE}, returns a dataframe of the inferred cell 
#'   type/lineage information instead.
#'
#' @importFrom Seurat AddMetaData as.SingleCellExperiment as.Seurat
#' @importFrom SingleR getReferenceDataset SingleR
#'
#' @export
#'
#' @examples
#' preds <- AssignCellType(pbmc_small, refset = "hpca", labels = "main_types",
#'   method = "single")
#'
#' preds2 <- AssignCellType(pbmc_small, refset = "hpca", labels = "types",
#'   method = "cluster", clusters = pbmc_small@meta.data$RNA_snn_res.1)
#'
AssignCellType <- function(scrna, refset = c("hpca", "blueprint_encode", 
  "immgen", "mouse.rnaseq"), labels = c("types", "main_types"), method = 
  c("single", "cluster"), clusters = NULL, assign = TRUE, n.cores = 1, ...) {

  # Arg matching.
  refset <- match.arg(refset)
  if (length(refset) > 1) {
    stop("No reference dataset specified.")
  }
  labels <- match.arg(labels)
  if (length(labels) > 1) {
    stop("'labels' must be specified.")
  }
  method <- match.arg(method)
  if (length(method) > 1) {
    stop("'method' must be specified.")
  }

  message("Downloading reference dataset: ", refset)
	dataset <- getReferenceDataset(dataset = refset)

  # Convert to SCE object and get common genes.
  sce <- as.SingleCellExperiment(scrna)
  common <- intersect(rownames(sce), rownames(dataset$data))
  sce <- sce[common,]
  sce <- sce[,colSums(counts(sce)) > 0]
  dataset$data <- dataset$data[common,]
  
  # Temporary workaround for sparse matrices throwing errors when 'cluster'
  # method used.
  if (method == "cluster") {
    counts(sce) <- as.matrix(counts(sce))
  }

  # Reference datasets are log2 transformed rather than natural log transformed
  # (Seurat). Scater normalization is log2. 
  sce <- scater::normalize(sce)

  message("Performing cell type inference.")
  annots <- SingleR(test = sce, training = dataset$data, 
    labels = hpca[[labels]], method = method, clusters = sce[[clusters]], 
    assay.type = 2, ...)

	if (isTRUE(assign)) {

	  message("Assigning inferred cell types.")

    if (method == "cluster") {
      annots$clusts <- rownames(annots)
      scrna[[sprintf("%s.%s.%s", clusters, refset, labels)]] <- 
        annots$labels[match(scrna@meta.data[[clusters]], 
          annots$clusts)]

    } else {
		  scrna[[sprintf("%s.%s", refset, labels)]] <- annots$labels
    }

	} else {
		scrna <- annots
	}

	return(scrna)
}

