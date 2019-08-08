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
#' If \code{outdir} is specified, the annotation results will be written to a
#' file named in \code{clusters.refset.labels.txt} format if 
#' \code{method="cluster"} or \code{refset.labels} format if 
#' \code{method="single"}.
#'
#' @param scrna Seurat object.
#' @param refset String indicating which reference dataset to download and use
#'   as the training set. 
#' @param labels If \code{TRUE}, use broad lineage-based labels 
#'   \code{main_types} for each reference sample. This is faster and can be more 
#'   informative depending on how specific your data is. If \code{FALSE}, more
#'   specific labels \code{types} will be used.
#' @param outdir Path to output directory for annotation scores, distributions,
#'   and heatmap.
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
#' @importFrom SingleCellExperiment counts
#' @importFrom utils write.table
#'
#' @export
#'
#' @examples
#' pbmc_small
#' preds <- AssignCellType(pbmc_small, refset = "hpca", labels = "main_types",
#'   method = "single")
#'
#' preds2 <- AssignCellType(pbmc_small, refset = "hpca", labels = "types",
#'   method = "cluster", clusters = "RNA_snn_res.1")
#'
AssignCellType <- function(scrna, refset = c("hpca", "blueprint_encode", 
  "immgen", "mouse.rnaseq"), labels = c("types", "main_types"), outdir = NULL,
  method = c("single", "cluster"), clusters = NULL, assign = TRUE, n.cores = 1, 
  ...) {

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
  dataset$data <- dataset$data[common,]

  # Reference datasets are log2 transformed rather than natural log transformed
  # (Seurat). Scater normalization is log2. 
  sce <- scater::logNormCounts(sce)

  message("Performing cell type inference.")
  annots <- SingleR(test = sce, training = dataset$data, 
    labels = dataset[[labels]], method = method, clusters = sce[[clusters]], 
    assay.type.test = "logcounts", assay.type.train = "logcounts", ...)

	if (isTRUE(assign)) {

	  message("Assigning inferred cell types.")
    label.dists <- 100 * (table(annots$labels) / 
      sum(table(annots$labels)))

    if (method == "cluster") {
      clusts <- rownames(annots)
      rownames(annots) <- NULL
      annots <- cbind(clusts, annots)

      if (!is.null(outdir)) {
        write.table(annots, file = sprintf("%s/%s.%s.%s.txt", outdir, clusters, 
          refset, labels), quote = FALSE, sep = "\t", row.names = FALSE)
        write.table(label.dists, file = sprintf("%s/%s.%s.%s.dist.txt", outdir, 
          clusters, refset, labels), quote = FALSE, sep = "\t", 
          row.names = FALSE)
        p <- plotScoreHeatmap(annots, clusters = annots$clusts, 
          silent = TRUE)
        pdf(sprintf("%s/%s.%s.%s.pdf", outdir, clusters, refset, labels))
        print(p)
        dev.off()
      }

      scrna[[sprintf("%s.%s.%s", clusters, refset, labels)]] <- 
        annots$labels[match(scrna@meta.data[[clusters]], 
          annots$clusters)]

    } else {
      cells <- rownames(annots)
      rownames(annots) <- NULL
      annots <- cbind(cells, annots)
      if (!is.null(outdir)) {
        write.table(annots, file = sprintf("%s/%s.%s.txt", outdir, refset, 
          labels), quote = FALSE, sep = "\t", row.names = FALSE)
        write.table(label.dists, file = sprintf("%s/%s.%s.dist.txt", outdir, 
          refset, labels), quote = FALSE, sep = "\t", row.names = FALSE)
        p <- plotScoreHeatmap(annots, clusters = sce[[clusters]], silent = TRUE)
        pdf(sprintf("%s/%s.%s.pdf", outdir, refset, labels))
        print(p)
        dev.off()
      }
		  scrna[[sprintf("%s.%s", refset, labels)]] <- annots$labels
    }

	} else {
		scrna <- annots
	}

	return(scrna)
}

