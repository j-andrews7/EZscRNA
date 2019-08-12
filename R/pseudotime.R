#' Infer principal root node of each partition
#'
#' \code{GetRootNodes} infers and returns the principal root node for each
#' partition. This will allow for the cells to be ordered according to the 
#' supplied variable (\code{order.var}) and start point \code{order.start}.
#'
#' @param cds \code{cell_data_set} object containing a trajectory graph, as from
#'   \code{\link{SeuratToCDS}}.
#' @param order.var String indicating \code{colData} column from which the 
#'   time points or progression order is defined.
#' @param order.start String indicating variable in \code{order.col} that
#'   represents the earliest time point or starting point.
#' @param partitions Character string indicating which partitions to calculate
#'   root nodes for. Done for all partitions if not provided.
#'
#' @return Character vector containing the names of the root nodes for each
#'   partition. If \code{partitions} parameter is set, only the root nodes for
#'   those partitions will be returned.
#'
#' @importFrom monocle3 principal_graph
#' @importFrom BiocGenerics which.max
#' @importMethodsFrom SummarizedExperiment colData
#'
#' @export
#'
GetRootNodes <- function(cds, order.var, order.start, partitions = NULL) {
  root.nodes <- c()

  if (is.null(partitions)) {
    partitions <- unique(cds$partitions)
  }

  for (i in partitions) {
    
    # Get cells in given partition.
    cell.ids <- which(colData(cds)[, "partitions"] == i & 
      colData(cds)[, order.var] == order.start)

    closest.vertex <-
      cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
    closest.vertex <- as.matrix(closest.vertex[colnames(cds), ])
    node <-
      igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
        (which.max(table(closest.vertex[cell.ids,]))))]

    root.nodes <- append(root.nodes, node)
  }

  return(root.nodes)
}


#' Convert a Seurat object to a cell data set object
#'
#' \code{SeuratToCDS} converts a \code{Seurat} object to a \code{cell_data_set}
#' object with a trajectory graph.
#'
#' @details
#' This function preserves the \code{UMAP} coordinates and clustering from the
#' \code{Seurat} object. This means any plots created from the 
#' \code{cell_data_set} will be comparable to those from 
#' \code{\link[Seurat]{DimPlot}}.
#'
#' @param scrna A Seurat object with UMAP mappings.
#' @param clusters String indicating column in \code{scrna@meta.data} that 
#'   contains cluster information to be utilized.
#' @return A \code{cell_data_set} object with trajectory information calculated.
#'   Retains the same UMAP mappings and cluster information as \code{scrna}.
#'
#' @importFrom monocle3 cluster_cells new_cell_data_set learn_graph partitions
#' @importMethodsFrom SingleCellExperiment reducedDims reducedDims<-
#' @importFrom SingleCellExperiment reducedDims
#' @import methods
#'
#' @export
#' 
SeuratToCDS <- function(scrna, clusters) {
  # Extract components to build CDS object.
  message("Creating cell_data_set object.")
  counts.data <- as(as.matrix(scrna@assays$RNA@data), 'sparseMatrix')
  cell_metadata <- scrna[[]]
  gene_annotation <- data.frame(gene_short_name = row.names(counts.data), 
    row.names = row.names(counts.data))

  # Build CDS object.
  cds <- new_cell_data_set(counts.data,
                           cell_metadata = cell_metadata,
                           gene_metadata = gene_annotation)

  # Assign UMAP coordinates to CDS object.
  reducedDims(cds)@listData[["UMAP"]] <- 
    scrna@reductions[["umap"]]@cell.embeddings

  # Cluster cells. Only done to get partitions that aren't provided by Seurat.
  cds <- cluster_cells(cds)

  # Map clusters to CDS object.
  list_cluster <- scrna@meta.data[[clusters]]
  names(list_cluster) <- scrna@assays[["RNA"]]@data@Dimnames[[2]]

  cds@clusters@listData[["UMAP"]][["clusters"]] <- list_cluster

  # Could be a place-holder, but essentially fills out louvain parameters.
  cds@clusters@listData[["UMAP"]][["louvain_res"]] <- "NA"

  # Assign UMAP coordinates to CDS object.
  reducedDims(cds)@listData[["UMAP"]] <- 
    scrna@reductions[["umap"]]@cell.embeddings

  message("Creating trajectory graph.")
  cds <- learn_graph(cds)

  cds$partitions <- partitions(cds)

  return(cds)
} 