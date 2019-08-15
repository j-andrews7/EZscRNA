#' Creates basic QC plots
#'
#' \code{RunQC} saves 3 QC plots showing gene counts, read counts, and percent 
#' mitochondrial reads per cell to help determine filters. It returns a
#' Seurat object with percent mitochondrial reads added to the \code{meta.data}.
#'
#' @param scrna Seurat object.
#' @param outdir Path to output directory for QC plots. Will not plot if not 
#'   set.
#' @return Seurat object with percent mitochondrial reads added to the
#'   \code{meta.data} for each cell.
#'
#' @importFrom Seurat PercentageFeatureSet VlnPlot
#'
#' @export
#'
#' @examples
#' library(Seurat)
#' pbmc_small <- RunQC(pbmc_small)
#'
RunQC <- function(scrna, outdir = NULL) {
  # Low-quality or dying cells often have mitochondrial contamination.
  scrna[["percent.mt"]] <- PercentageFeatureSet(scrna, pattern = "^MT-")

  if (!is.null(outdir)) {
    message("Creating QC plots.")
    pdf(sprintf("%s/Vln.QC.Metrics.pdf", outdir), height = 12, width = 8, 
    	useDingbats = FALSE)
    p <- VlnPlot(scrna, features = c("nFeature_RNA", "nCount_RNA", 
  		"percent.mt"), ncol = 1)
    print(p)
    dev.off()
  }

  return(scrna)
}


#' Normalize counts and score cell cycle for each cell
#' 
#' \code{NormScoreCC} returns a Seurat object with normalized counts and adds 
#' cell cycle scores for each gene based on Seurat's cell cycle gene lists.
#'
#' The Seurat authors state (https://github.com/satijalab/seurat/issues/1679) 
#' that counts should always be normalized before cell cycle or module scoring.
#' This is particularly important if one is using the \code{SCTranform} function
#' for data normalization, scaling, and therefore, regression. 
#'
#' @param scrna Seurat object to score cell cycle genes for each cell.
#' @param skip.sct Boolean indicating whether to skip \code{SCTransform} call.
#'   Useful for integrated objects.
#' @return Seurat object with cell cycle scores ('S.Score', 'G2M.Score') and 
#'   'Phase' added to \code{meta.data} for each cell.
#'
#' @importFrom Seurat NormalizeData CellCycleScoring
#'
#' @export
#'
#' @examples
#' \dontrun{
#' library(Seurat)
#' pbmc_small <- NormScoreCC(pbmc_small)
#' }
#'
NormScoreCC <- function(scrna, skip.sct = NULL) {
	message("Scoring cell cycle genes.")
	# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with 
	# Seurat. We can segregate this list into markers of G2/M and of S phase.
	s.genes <- cc.genes$s.genes
	g2m.genes <- cc.genes$g2m.genes

  if (!isTRUE(skip.sct)) {
	  scrna <- SCTransform(scrna)
  }

	cc.seurat <- CellCycleScoring(scrna, s.features = s.genes, 
		g2m.features = g2m.genes)
	message("Cell cycle scoring complete.")

	return(cc.seurat)
}


#' Exploratory data analysis plots
#'
#' \code{BatchCCEDA} creates a number of plots to determine the number of PCs
#' to use for PCA/clustering and whether or not cell cycle scores and batch
#' effects should be addressed. Runs and plots an ElbowPlot to determine PCs
#' for later use. Runs and plots PCA for cell cycle genes to show their impact.
#'
#' @details
#' Supplying \code{vars} will plot a PCA from the variable genes for each 
#' variable. It will also calculate and create a density plot of the variance 
#' explained by each variable across all genes using 
#' \code{\link[scater]{plotExplanatoryVariables}}.
#'
#' @param scrna Seurat object.
#' @param outdir Path to output directory for plots.
#' @param npcs Number of PCs to use for PCA and ElbowPlot. 
#' @param vars Character vector indicating \code{meta.data} columns to be 
#'   investigated for batch effects and variance contributions.
#' @param skip.sct Boolean indicating whether to skip \code{SCTransform} call.
#'   Useful for integrated objects.
#' @return A Seurat object with a PCA for cell cycle genes stored with
#'   \code{reduction.name = "cc"}. 
#'
#' @importFrom Seurat SCTransform RunPCA ElbowPlot DimPlot 
#'   as.SingleCellExperiment
#' @importFrom grDevices dev.off pdf
#' @importFrom scater plotExplanatoryVariables
#'
#' @export
#'
#' @examples
#' library(Seurat)
#' \dontrun{
#' pbmc_small <- RunQC(pbmc_small)
#' pbmc_small <- NormScoreCC(pbmc_small)
#'
#' # Skip SCT normalization - should only be done if object has been integrated.
#' pbmc_small <- BatchCCEDA(pbmc_small, skip.sct = TRUE)
#'
#' # Can explore other variables as well.
#' pbmc_small <- BatchCCEDA(pbmc_small, skip.sct = TRUE, vars = "group")
#' }
#'
BatchCCEDA <- function(scrna, outdir = ".", npcs = 50, vars = NULL, 
  skip.sct = NULL) {

  # A list of cell cycle markers, from Tirosh et al, 2015, is loaded with 
  # Seurat.  We can segregate this list into markers of G2/M and S phase.
  s.genes <- cc.genes$s.genes
  g2m.genes <- cc.genes$g2m.genes

  # Scale and return all genes so cell cycle gene PCAs won't lie to use.
  if (!isTRUE(skip.sct)) {
    scrna <- SCTransform(scrna, vars.to.regress = c("percent.mt", "nCount_RNA"),
  	  return.only.var.genes = FALSE)
  }

  scrna <- RunPCA(scrna, npcs = npcs)

  # Elbowplot to determine number of PCs to use later.
  message("Creating ElbowPlot.")
  pdf(sprintf("%s/ElbowPlot.pdf", outdir), useDingbats = FALSE)
  p <- ElbowPlot(scrna, ndims = npcs)
  print(p)
  dev.off()

  # PCA for variable genes across vars.
  message("Creating PCA/variance plots using variable genes.")
  if(!is.null(vars)) {
    for (i in vars) {
  		pdf(sprintf("%s/PCA.%s.NoRegression.pdf", outdir, i), height = 5, 
        width = 7, useDingbats = FALSE)
  		p <- DimPlot(scrna, group.by = i, pt.size = 0.3)
  		print(p)
  		dev.off()
    }
    sce <- as.SingleCellExperiment(scrna)
    var.feats <- VariableFeatures(scrna)
    p1 <- plotExplanatoryVariables(sce, variables = vars) + ggtitle("All Genes")
    p2 <- plotExplanatoryVariables(sce, variables = vars, 
      subset_row = var.feats) + 
      ggtitle(paste0("Top ", length(var.feats), " Variable Genes"))

    pdf(sprintf("%s/DensityPlot.variance.pdf", outdir), height = 5, 
        width = 7, useDingbats = FALSE) 
    print(p1)
    print(p2)
    dev.off()
  }

  # Now take a look at the PCA for the cell cycle genes.
  message("Performing cell cycle PCA.")
  scrna <- RunPCA(scrna, npcs = npcs, features = c(s.genes, g2m.genes),
  	reduction.name = "cc")
  pdf(sprintf("%s/PCA.CellCycle.NoRegression.pdf", outdir), height = 5, 
  	width = 7, useDingbats = FALSE)
  p <- DimPlot(scrna, group.by = "Phase", pt.size = 0.3, reduction = "cc")
  print(p)
  dev.off()

  return(scrna)
}