#' Creates basic QC plots
#'
#' \code{RunQC} saves 3 QC plots showing gene counts, read counts, and percent 
#' mitochondrial reads per cell to help determine filters. It returns a
#' Seurat object with percent mitochondrial reads added to the \code{meta.data}.
#'
#' @param scrna Seurat object.
#' @param outdir Path to output directory.
#' @return Seurat object with percent mitochondrial reads added to the
#'   \code{meta.data} for each cell.
#'
#' @importFrom Seurat PercentageFeatureSet VlnPlot
#'
#' @export
#'
RunQC <- function(scrna, outdir = ".") {
  # Low-quality or dying cells often have mitochondrial contamination.
  scrna[["percent.mt"]] <- PercentageFeatureSet(scrna, pattern = "^MT-")

  message("Creating QC plots.")
  pdf(sprintf("%s/Vln.QC.Metrics.pdf", outdir), height = 12, width = 8, 
  	useDingbats = FALSE)
  p <- VlnPlot(scrna, features = c("nFeature_RNA", "nCount_RNA", 
		"percent.mt"), ncol = 1)
  print(p)
  dev.off()

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
#' @return Seurat object with cell cycle scores ('S.Score', 'G2M.Score') and 
#'   'Phase' added to \code{meta.data} for each cell.
#'
#' @importFrom Seurat NormalizeData CellCycleScoring
#'
#' @export
#'
NormScoreCC <- function(scrna) {
	message("Scoring cell cycle genes.")
	# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with 
	# Seurat. We can segregate this list into markers of G2/M and of S phase.
	s.genes <- cc.genes$s.genes
	g2m.genes <- cc.genes$g2m.genes

	scrna <- NormalizeData(scrna)

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
#' PCA on variable features can also be plotted by batch to view potential batch
#' effects.
#'
#' @param scrna Seurat object.
#' @param outdir Path to output directory for plots.
#' @param npcs Number of PCs to use for PCA and ElbowPlot. 50 by default.
#' @param batch String indicating \code{meta.data} column to be investigated for
#'   batch effects. 
#' @return A Seurat object with a PCA for cell cycle genes stored with
#'   \code{reduction.name = "cc"}. 
#'
#' @importFrom Seurat SCTransform RunPCA ElbowPlot DimPlot
#' @importFrom grDevices dev.off pdf
#'
#' @export
#'
BatchCCEDA <- function(scrna, outdir, npcs = 50, batch = NULL) {

  # A list of cell cycle markers, from Tirosh et al, 2015, is loaded with 
  # Seurat.  We can segregate this list into markers of G2/M and S phase.
  s.genes <- cc.genes$s.genes
  g2m.genes <- cc.genes$g2m.genes

  # Scale and return all genes so cell cycle gene PCAs won't lie to use.
  scrna <- SCTransform(scrna, vars.to.regress = "percent.mt",
  	return.only.var.genes = FALSE)
  scrna <- RunPCA(scrna, npcs = npcs)

  # Elbowplot to determine number of PCs to use later.
  pdf(sprintf("%s/ElbowPlot.pdf", outdir), useDingbats = FALSE)
  p <- ElbowPlot(scrna, ndims = npcs)
  print(p)
  dev.off()

  # Take a look at PCA for variable genes across batch.
  if(!is.null(batch)) {
		pdf(sprintf("%s/PCA.Batch.NoRegression.pdf", outdir), height = 5, width = 7,
			useDingbats = FALSE)
		p <- DimPlot(scrna, group.by = batch, pt.size = 0.3)
		print(p)
		dev.off()
  }

  # Now take a look at the PCA for the cell cycle genes.
  scrna <- RunPCA(scrna, npcs = npcs, features = c(s.genes, g2m.genes),
  	reduction.name = "cc")
  pdf(sprintf("%s/PCA.CellCycle.NoRegression.pdf", outdir), height = 5, 
  	width = 7, useDingbats = FALSE)
  p <- DimPlot(scrna, group.by = "Phase", pt.size = 0.3, reduction = "cc")
  print(p)
  dev.off()

  return(scrna)
}