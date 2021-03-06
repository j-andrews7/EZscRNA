% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/eda.R
\name{BatchCCEDA}
\alias{BatchCCEDA}
\title{Exploratory data analysis plots}
\usage{
BatchCCEDA(scrna, outdir = ".", npcs = 50, vars = NULL,
  skip.sct = NULL)
}
\arguments{
\item{scrna}{\linkS4class{Seurat} object.}

\item{outdir}{Path to output directory for plots.}

\item{npcs}{Number of PCs to use for PCA and ElbowPlot.}

\item{vars}{Character vector indicating \code{meta.data} columns to be 
investigated for batch effects and variance contributions.}

\item{skip.sct}{Boolean indicating whether to skip \code{SCTransform} call.
Useful for integrated objects.}
}
\value{
A \linkS4class{Seurat} object with a PCA for cell cycle genes stored 
  with \code{reduction.name = "cc"}.
}
\description{
\code{BatchCCEDA} creates a number of plots to determine the number of PCs
to use for PCA/clustering and whether or not cell cycle scores and batch
effects should be addressed. Runs and plots an ElbowPlot to determine PCs
for later use. Runs and plots PCA for cell cycle genes to show their impact.
}
\details{
Supplying \code{vars} will plot a PCA from the variable genes for each 
variable. It will also calculate and create a density plot of the variance 
explained by each variable across all genes using 
\code{\link[scater]{plotExplanatoryVariables}}.
}
\examples{
library(Seurat)
\dontrun{
pbmc_small <- RunQC(pbmc_small)
pbmc_small <- NormScoreCC(pbmc_small)

# Skip SCT normalization - should only be done if object has been integrated.
pbmc_small <- BatchCCEDA(pbmc_small, skip.sct = TRUE)

# Can explore other variables as well.
pbmc_small <- BatchCCEDA(pbmc_small, skip.sct = TRUE, vars = "group")
}

}
\author{
Jared Andrews
}
