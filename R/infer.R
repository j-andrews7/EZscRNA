#' Infers and assigns nearest cell type for each cell using SingleR
#'
#' \code{AssignCellType} performs correlation-based cell type inference using a 
#' reference dataset via \code{\link[SingleR]{SingleR}}. 
#'
#' @details
#' Reference datasets available for use by this function via 
#' \code{\link[SingleR]{getReferenceDataset}} include:
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
#' \code{method="single"}. Additionally, a heatmap will be made from the
#' annotation results via \code{\link[SingleR]{plotScoreHeatmap}}.
#'
#' @param scrna Seurat object.
#' @param refset Reference dataset as retrieved by 
#'   \code{\link[SingleR]{SingleR}}.
#' @param labels String indicating whether to use broad lineage-based labels 
#'   \code{labels="main_types"} for each reference sample or more fine-grained 
#'   labels (\code{labels="types"}. The former is quicker and can be informative
#'   enough if your sample has many cell types. The latter is best-suited for 
#'   purified cell types or if particular cellular subtypes are important.
#' @param outdir Path to output directory for annotation scores, distributions,
#'   and heatmap.
#' @param method String specifying whether annotation should be applied to 
#'   each single cell \code{method = "single"} or aggregated into cluster-level 
#'   profiles \code{method = "cluster"} prior to annotation.
#' @param clusters String defining the \code{meta.data} column in \code{scrna}
#'   that specify cluster identities for each cell. Only required if 
#'   \code{method = "cluster"}. If provided with \code{method = "single"}, will
#'   be used as an additional label in heatmap plotting.
#' @param integrated Boolean indicating whether Seurat object was integrated via
#'   \code{\link{SimpleIntegration}}. This must be set to \code{TRUE} if so,
#'   as the normalized counts are not stored in the "integrated" assay.
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
#' @importFrom Seurat AddMetaData as.SingleCellExperiment
#' @importFrom SingleR SingleR plotScoreHeatmap
#' @importFrom SingleCellExperiment counts
#' @importFrom utils write.table
#'
#' @export
#'
#' @examples
#' library("Seurat")
#' library("SingleR")
#' pbmc <- pbmc_small
#' hpca <- getReferenceDataset(dataset = "hpca")
#' preds <- AssignCellType(pbmc, refset = hpca, labels = "main_types",
#'   method = "single")
#'
#' preds2 <- AssignCellType(pbmc, refset = hpca, labels = "types",
#'   method = "cluster", clusters = "RNA_snn_res.1")
#'
AssignCellType <- function(scrna, refset, labels = c("types", "main_types"), 
  outdir = NULL, method = c("single", "cluster"), clusters = NULL, 
  integrated = FALSE, assign = TRUE, n.cores = 1, ...) {

  # Arg matching.
  labels <- match.arg(labels)
  if (length(labels) > 1) {
    stop("'labels' must be specified.")
  }
  method <- match.arg(method)
  if (length(method) > 1) {
    stop("'method' must be specified.")
  }

  # Convert to SCE object and get common genes.
  if (isTRUE(integrated)) {
    sce <- as.SingleCellExperiment(scrna, assay = "SCT")
  } else {
    sce <- as.SingleCellExperiment(scrna)
  }

  common <- intersect(rownames(sce), rownames(refset$data))
  sce <- sce[common,]
  refset$data <- refset$data[common,]

  # Get clusters so subsetting by NULL isn't an issue.
  if (!is.null(clusters)) {
    clusts <- sce[[clusters]]
  } else {
    clusts <- NULL
  }

  # Reference datasets are log2 transformed rather than natural log transformed
  # (Seurat). Scater normalization is log2. 
  sce <- scater::logNormCounts(sce)

  message("Performing cell type inference.")
  annots <- SingleR(test = sce, ref = refset$data, 
    labels = refset[[labels]], method = method, clusters = clusts,
    ...)

	if (isTRUE(assign)) {

	  message("Assigning inferred cell types.")
    label.dists <- 100 * (table(annots$labels) / 
      sum(table(annots$labels)))

    if (method == "cluster") {
      clusts <- rownames(annots)
      rownames(annots) <- NULL
      annots <- cbind(clusts, annots)

      if (!is.null(outdir)) {
        # Create output tables/plots.
        write.table(annots, file = sprintf("%s/%s.%s.%s.txt", outdir, clusters, 
          refset$name, labels), quote = FALSE, sep = "\t", row.names = FALSE)
        write.table(label.dists, file = sprintf("%s/%s.%s.%s.dist.txt", outdir, 
          clusters, refset$name, labels), quote = FALSE, sep = "\t", 
          row.names = FALSE)
        p <- plotScoreHeatmap(annots, clusters = annots$clusts, 
          silent = TRUE)
        pdf(sprintf("%s/%s.%s.%s.pdf", outdir, clusters, refset$name, labels))
        print(p)
        dev.off()
      }

      # Add labels as column to meta.data.
      scrna[[sprintf("%s.%s.%s", clusters, refset$name, labels)]] <- 
        annots$labels[match(scrna@meta.data[[clusters]], annots$clusts)]

    } else {
      cells <- rownames(annots)
      rownames(annots) <- NULL
      annots <- cbind(cells, annots)

      if (!is.null(outdir)) {
        if(!is.null(clusters)) {
          clusts <- sce[[clusters]]
          clust.label <- clusters
        } else{
          clusts <- NULL
          clust.label <- "NoClusters"
        }
        # Create output tables/plots.
        write.table(annots, file = sprintf("%s/%s.%s.txt", outdir, refset$name, 
          labels), quote = FALSE, sep = "\t", row.names = FALSE)
        write.table(label.dists, file = sprintf("%s/%s.%s.dist.txt", outdir, 
          refset$name, labels), quote = FALSE, sep = "\t", row.names = FALSE)
        if (nrow(annots) < 65500) {
          p <- plotScoreHeatmap(annots, clusters = clusts, silent = TRUE)
          pdf(sprintf("%s/%s.%s.%s.%s.pdf", outdir, method, refset$name, labels,
            clust.label))
          print(p)
          dev.off()
        } else {
          message(paste0("Skipping `plotScoreHeatmap()`, as it can't cluster", 
            " more than 65000 cells."))
        }
      }
		  scrna[[sprintf("%s.%s", refset$name, labels)]] <- annots$labels
    }

	} else {
		scrna <- annots
	}

	return(scrna)
}

