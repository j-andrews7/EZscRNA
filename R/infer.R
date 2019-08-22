#' Infers and assigns nearest cell type for each cell using SingleR
#'
#' \code{AssignCellType} performs correlation-based cell type inference using a 
#' reference dataset via \code{\link[SingleR]{SingleR}}. 
#'
#' @details
#' Reference datasets available for use by this function include:
#' \describe{
#'   \item{"HPCA"}{Human Primary Cell Atlas (HPCA): A collection of Gene 
#'     Expression Omnibus (GEO datasets), which contains 713 human microarray 
#'     samples classified to 38 main cell types and further annotated to 169
#'     subtypes. Download from ExperimentHub as a 
#'     \linkS4class{SummarizedExperiment} object via 
#'     \code{\link[SingleR]{HumanPrimaryCellAtlasData}}. See dataset details via 
#'     \code{?HumanPrimaryCellAtlasData}.}
#'   \item{"Blueprint_Encode"}{Blueprint + ENCODE datasets: Blueprint 
#'     Epigenomics, 144 RNA-seq pure human immune samples annotated to 28 cell 
#'     types. ENCODE, 115 RNA-seq pure stroma and immune samples annotated 
#'     to 17 cell types. Altogether, 259 samples with 43 cell types. Download 
#'     from ExperimentHub as a \linkS4class{SummarizedExperiment} object via 
#'     \code{\link[SingleR]{BlueprintEncodeData}}. See dataset details via 
#'     \code{?BlueprintEncodeData}.}
#'   \item{"ImmGen"}{Immunological Genome Project (ImmGen): 830 microarray 
#'     samples, which we classified to 20 main cell types and further
#'     annotated to 253 subtypes. Download from ExperimentHub as a 
#'     \linkS4class{SummarizedExperiment} object via 
#'     \code{\link[SingleR]{ImmGenData}}. See dataset details via 
#'     \code{?ImmGeneData}.}
#'   \item{"MouseRNAseq"}{A dataset of 358 mouse RNA-seq samples annotated to
#'     28 cell types. This dataset was collected, processed, and shared 
#'     courtesy of Bérénice Benayoun. This data set is especially useful for 
#'     brain-related samples. Download from ExperimentHub as a 
#'     \linkS4class{SummarizedExperiment} object via 
#'     \code{\link[SingleR]{MouseRNAseqData}}. See dataset details via 
#'     \code{?MouseRNAseqData}.}
#' }
#'
#' If \code{outdir} is specified, the annotation results will be written to a
#' file named in \code{clusters.refset.labels.txt} format if 
#' \code{method="cluster"} or \code{refset.labels} format if 
#' \code{method="single"}. Additionally, a heatmap will be made from the
#' annotation results via \code{\link[SingleR]{plotScoreHeatmap}} if there are
#' fewer than 65500 cells (the hclust method used fails with more cells than 
#' that). 
#'
#' @param scrna Seurat object.
#' @param refsets List of reference dataset(s). Multiple can be provided.
#' @param labels String indicating whether to use broad lineage-based labels 
#'   \code{labels="label.main"} for each reference sample or more fine-grained 
#'   labels (\code{labels="label.fine"}. The former is quicker and can be 
#'   informative enough if your sample has many cell types. The latter is 
#'   best-suited for purified cell types or if particular cellular subtypes are 
#'   important.
#' @param outdir Path to output directory for annotation scores, distributions,
#'   and heatmap.
#' @param method String or character vector specifying whether annotation should 
#'   be applied to each single cell \code{method = "single"} or aggregated into 
#'   cluster-level profiles \code{method = "cluster"} prior to annotation. If
#'   both are supplied, predictions will be performed for both methods.
#' @param clusters String or character vector defining the \code{meta.data} 
#'   column(s) in \code{scrna} that specify cluster identities for each cell. 
#'   Only required if \code{method = "cluster"}. A character vector can be 
#'   provided if multiple annotations with different clusters is wanted.
#'   
#'   If provided with \code{method = "single"}, will be used as an additional 
#'   label in heatmap plotting. 
#' @param singler.params List of arguments to be passed to 
#'   \code{\link[SingleR]{SingleR}}.
#' @return A \code{Seurat} object with inferred
#'   cell type information in the \code{meta.data} slot named in 
#'   \code{refset.labels} format. If \code{method = "cluster"}, the resulting
#'   \code{meta.data} column will be named in \code{clusters.refset.labels} 
#'   format. 
#'
#' @importFrom utils write.table
#'
#' @export
#'
#' @examples
#' library("Seurat")
#' library("SingleR")
#' pbmc <- pbmc_small
#' # Download reference data from ExperimentHub.
#' hpca <- HumanPrimaryCellAtlasData()
#' # Always give the reference set a name.
#' metadata(hpca)$name <- "HPCA"
#' pbmc.pred <- AssignCellType(pbmc, refsets = list(hpca), 
#'   labels = "label.main", method = "single")
#'
#' pbmc.clusters.pred <- AssignCellType(pbmc, refsets = list(hpca), 
#'   labels = "label.fine", method = "cluster", clusters = "RNA_snn_res.1")
#'
#' \dontrun{
#' # Use both broad and fine labels, 'Single' and 'cluster' method, multiple 
#' # reference sets, and multiple cluster annotations.
#' bp <- BlueprintEncodeData()
#' metadata(bp)$name <- "Blueprint_Encode"
#' pbmc.all.anno <- AssignCellType(pbmc, refsets = list(hpca, bp), 
#'   clusters = c("RNA_snn_res.0.8", "RNA_snn_res.1"))
#' }
#'
#' @author Jared Andrews
#'
AssignCellType <- function(scrna, refsets, labels = c("label.main", 
  "label.fine"), outdir = NULL, method = c("single", "cluster"), 
  clusters = NULL, singler.params = NULL) {

  # Arg matching & checking.
  labels <- match.arg(labels, several.ok = TRUE)
  method <- match.arg(method, several.ok = TRUE)

  for (ref in refsets) {
    if(is.null(ref@metadata$name)) {
      stop("The refset object must have a name", 
        " (e.g. hpca@metadata$name <- 'HPCA'")
    }
  }

  if ((is.null(clusters)) && (method == "cluster")) {
    stop("'clusters' cannot be NULL when 'method=\"cluster\"'.")
  }

  # Convert to SingleCellExperiment.
  sce <- Seurat::as.SingleCellExperiment(scrna, assay = "RNA")

  # Reference datasets are log2 transformed rather than natural log transformed
  # (Seurat). Scater normalization is log2. 
  sce <- scater::logNormCounts(sce)

  for (ref in refsets) {
    message("Get common genes in ", ref@metadata$name)
    common <- intersect(rownames(sce), rownames(ref))
    sce.com <- sce[common,]
    ref <- ref[common,]

    # This is admittedly a touch gross.
    for (meth in method) {
      for (lab in labels) {
        if (!is.null(clusters)) {
          for (c in clusters) {
            scrna <- PerformCellInference(scrna = scrna, sce = sce.com, 
              method = meth, refset = ref, labels = lab, 
              clusts.name = c, clusts = sce[[c]], outdir = outdir, 
              singler.params = singler.params)
          }
        } else {
          scrna <- PerformCellInference(scrna = scrna, sce = sce.com, 
            method = meth, refset = ref, labels = lab, outdir = outdir,
            singler.params = singler.params)
        }
      }
    }
  }

	return(scrna)
}

#' Performs SingleR inference
#'
#' \code{PerformCellInference} performs correlation-based cell type inference 
#' using a reference dataset via \code{\link[SingleR]{SingleR}}. 
#'
#' @param scrna \code{Seurat} object.
#' @param sce \code{SingleCellExperiment} object.
#' @param method String or character vector specifying whether annotation should 
#'   be applied to each single cell \code{method = "single"} or aggregated into 
#'   cluster-level profiles \code{method = "cluster"} prior to annotation. If
#'   both are supplied, predictions will be performed for both methods.
#' @param refset Reference dataset as \linkS4class{SummarizedExperiment}. 
#' @param labels String indicating whether to use broad lineage-based labels 
#'   \code{labels="label.main"} for each reference sample or more fine-grained 
#'   labels (\code{labels="label.fine"}. 
#' @param outdir Path to output directory for annotation scores, distributions,
#'   and heatmap.
#' @param clusts.name String indicating clusters column.
#' @param clusts Vector of cluster assignments for each cell.
#' @param singler.params List of arguments to be passed to 
#'   \code{\link[SingleR]{SingleR}}.
#' @return A \code{Seurat} object with inferred
#'   cell type information in the \code{meta.data} slot named in 
#'   \code{refset.labels} format. If \code{method = "cluster"}, the resulting
#'   \code{meta.data} column will be named in \code{clusters.refset.labels} 
#'   format. 
#'
#' @importFrom utils write.table
#' @importFrom grDevices dev.off pdf
#' @importFrom S4Vectors metadata
#'
#' @author Jared Andrews
#'
PerformCellInference <- function(scrna, sce, method, refset, labels, 
  outdir, clusts.name = "NoClusters", clusts = NULL, singler.params = NULL) {

  message("Performing cell type inference.")
  
  start <- Sys.time()
  annots <- do.call(SingleR::SingleR, c(list(test = sce, ref = refset, 
    labels = refset[[labels]], method = method, clusters = clusts), 
    singler.params))
  end <- Sys.time()

  elapsed <- end - start
  message("Took: ", round(elapsed, digits = 3), " seconds with:\n", 
    paste0("labels: ", labels, "; method: ", method, "; ref: ", 
      metadata(refset)$name, "; clusters: ", clusts.name))

  label.dists <- 100 * (table(annots$labels) / 
    sum(table(annots$labels)))

  if (method == "cluster") {
    clusters <- rownames(annots)
    rownames(annots) <- NULL
    annots <- cbind(clusters, annots)

    if (!is.null(outdir)) {
      # Create output tables/plots.
      write.table(annots, file = sprintf("%s/%s.%s.%s.txt", outdir, clusts.name, 
        metadata(refset)$name, labels), quote = FALSE, sep = "\t", 
        row.names = FALSE)
      write.table(label.dists, file = sprintf("%s/%s.%s.%s.dist.txt", outdir, 
        clusts.name, metadata(refset)$name, labels), quote = FALSE, sep = "\t", 
        row.names = FALSE)
      p <- SingleR::plotScoreHeatmap(annots, clusters = annots$clusts, 
        silent = TRUE)
      pdf(sprintf("%s/%s.%s.%s.pdf", outdir, clusts.name, metadata(refset)$name, 
        labels))
      print(p)
      dev.off()
    }

    # Add labels as column to meta.data.
    scrna[[sprintf("%s.%s.%s", clusts.name, metadata(refset)$name, labels)]] <- 
      annots$labels[match(scrna@meta.data[[clusts.name]], annots$clusters)]

  } else {
    cells <- rownames(annots)
    rownames(annots) <- NULL
    annots <- cbind(cells, annots)

    if (!is.null(outdir)) {

      # Create output tables/plots.
      write.table(annots, file = sprintf("%s/%s.%s.txt", outdir, 
        metadata(refset)$name, labels), quote = FALSE, sep = "\t", 
        row.names = FALSE)
      write.table(label.dists, file = sprintf("%s/%s.%s.dist.txt", outdir, 
        metadata(refset)$name, labels), quote = FALSE, sep = "\t", 
        row.names = FALSE)

      # pheatmap clustering can't handle tons of cells.
      if (nrow(annots) < 65500) {
        pdf(sprintf("%s/%s.%s.%s.%s.pdf", outdir, method, refset@metadata$name, 
          labels, clusts.name))
        p <- SingleR::plotScoreHeatmap(annots, clusters = clusts, 
          silent = TRUE)
        print(p)
        dev.off()
      } else {
        message(paste0("Skipping `plotScoreHeatmap()`, as it can't cluster", 
          " more than 65000 cells."))
      }
    }
    scrna[[sprintf("%s.%s", metadata(refset)$name, labels)]] <- annots$labels
  }

  return(scrna)
}
