#' Infers and assigns nearest cell type for each cell using SingleR
#'
#' \code{InferCellType} performs correlation-based cell type inference using a 
#' reference dataset via \code{\link[SingleR]{SingleR}}. 
#'
#' @details
#' Reference datasets available for use by this function include those provided
#' by SingleR:
#' \itemize{
#'   \item \code{\link[SingleR]{HumanPrimaryCellAtlasData}}
#'   \item \code{\link[SingleR]{BlueprintEncodeData}}
#'   \item \code{\link[SingleR]{DatabaseImmuneCellExpressionData}}
#'   \item \code{\link[SingleR]{NovershternHematopoieticData}}
#'   \item \code{\link[SingleR]{MonacoImmuneData}}
#'   \item \code{\link[SingleR]{ImmGenData}}
#'   \item \code{\link[SingleR]{MouseRNAseqData}}
#' }
#'
#' If \code{outdir} is specified, the annotation results will be written to a
#' file named in \code{clusters.refset.labels.txt} format if 
#' \code{method="cluster"} or \code{refset.labels.txt} format if 
#' \code{method="single"}. Label distributions will be written to files named in
#' \code{refset.labels.dist.txt} and \code{refset.labels.pruned.dist.txt} 
#' format. Additionally, a heatmap will be made from the annotation results via 
#' \code{\link[SingleR]{plotScoreHeatmap}} if there are fewer than 65500 cells 
#' (the hclust method used fails with more cells than that). 
#'
#' @param scrna \linkS4class{Seurat} object.
#' @param refs List of reference dataset(s). Multiple can be provided.
#' @param labels List of vectors containing vectors for labels for each sample
#'   in each reference dataset.
#' @param outdir Path to output directory for results table(s).
#' @param clusters String or character vector defining the \code{meta.data} 
#'   column(s) in \code{scrna} that specify cluster identities for each cell. 
#'   A character vector can be provided if multiple annotations with different 
#'   clusters is wanted.
#'   
#'   If provided with \code{method = "single"}, clusters will be used as an 
#'   additional label in heatmap plotting. 
#' @param ... Extra parameters to be passed to 
#'   \code{\link[SingleR]{SingleR}}.
#' @return A \linkS4class{Seurat} object with inferred
#'   cell type information in the \code{meta.data} slot named in 
#'   \code{refset.labels} format. If \code{clusters} is provided, the resulting
#'   \code{meta.data} column will be named in \code{clusters.refset.labels} 
#'   format. Pruned labels will also be added with \code{.pruned} prepended to
#'   the column name.
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
#' pbmc.pred <- InferCellType(pbmc, refs = list(HPCA=hpca), 
#'   labels = list(hpca$label.fine), method = "single")
#'
#' pbmc.clusters.pred <- InferCellType(pbmc, refsets = list(HPCA=hpca), 
#'   labels = list(hpca$label.fine), method = "cluster", clusters = "RNA_snn_res.1")
#'
#' @seealso \code{\link[SingleR]{SingleR}} for additional options.
#'
#' @author Jared Andrews
#'
InferCellType <- function(scrna, refs, labels, outdir = NULL, clusters = NULL, 
    recomp = TRUE, ...) {

    dir.create(file.path(outdir), recursive = TRUE, showWarnings = FALSE)
    
    sce <- Seurat::as.SingleCellExperiment(scrna, assay = "RNA")
    sce <- scater::logNormCounts(sce)
    
    if (recomp) {
        out.suf <- "recomp"
    } else {
        # Used to differentiate between common combining and recomputed combining.
        out.suf <- "common"
    }

    pred <- SingleR::SingleR(test = sce, ref = refs, labels = labels, method = "single", 
        recompute = recomp, ...)

    pred$cells <- rownames(pred)

    if (!is.null(outdir)) {
        write.table(pred, file = sprintf("%s/CellPredictions.%s.txt", outdir, out.suf), 
            sep = "\t", quote = FALSE, row.names = FALSE)
    }

    scrna[[sprintf("%s.pruned.labels", out.suf)]] <- pred$pruned.labels
    scrna[[sprintf("%s.labels", out.suf)]] <- pred$labels
    
    if (!is.null(clusters)) {
        for (clust in clusters) {

            pred <- SingleR::SingleR(test = sce, ref = refs, labels = labels, 
                method = "cluster", clusters = sce[[clust]], recompute = recomp, ...)

            pred$clusters <- rownames(pred)
            scrna[[sprintf("%s.%s.pruned.pred",clust, out.suf)]] <- pred$pruned.labels[match(scrna[[]][[clust]], rownames(pred))]
            scrna[[sprintf("%s.%s.pred",clust, out.suf)]] <- pred$labels[match(scrna[[]][[clust]], rownames(pred))]

            if (!is.null(outdir)) {
                write.table(pred, file = sprintf("%s/%s.%s.ClusterPredictions.txt", outdir, clust, out.suf), 
                    sep = "\t", quote = FALSE, row.names = FALSE)
            }
        }
    }
    
    return(scrna)
}
