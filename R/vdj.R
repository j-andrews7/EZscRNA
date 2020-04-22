#' Add 10X clonotype data to seurat object
#'
#' \code{AddClonotype} adds clonotype info from matched 10X VDJ sequencing to 
#' the metadata of a given \linkS4class{Seurat} object.
#'
#' This function is admittedly rough and will be rewritten in the future. It
#' does not include specific VDJ genes for each cell, rather just using the
#' final amino acid sequence for inter-sample comparison.
#'
#' @param vdj data.frame containing contig info.
#' @param scrna \linkS4class{Seurat} object.
#' @param freq.cutoff Clonotype frequency cutoff for retaining cdr3 sequence
#'   in \code{cdr3_group} column rather than grouping it as "Other".
#' @return \linkS4class{Seurat} object with clonotype data added to the 
#'   metadata for each cell.
#'
#' @importFrom Seurat AddMetaData
#' @importFrom utils read.csv
#' @importFrom CellaRepertorium ContigCellDB_10XVDJ canonicalize_cell
#'
#' @export
#'
AddClonotype <- function(vdj, scrna, freq.cutoff = 0.01) {
  
  vdj$barcode <- gsub("-1", "", vdj$barcode)
  vdj <- ContigCellDB_10XVDJ(vdj)
    
  can.vdj <- canonicalize_cell(vdj, contig_fields = c('umis', 'reads', 
    'chain', 'v_gene', 'd_gene', 'j_gene', 'cdr3'))
    
  # Reorder so barcodes are first column and set them as rownames.
  vdj <- data.frame(can.vdj$cell_tbl)
  rownames(vdj) <- vdj[, 8]
  vdj[, 8] <- NULL

  # Add to the Seurat object's metadata and convert to string.
  clono.seurat <- AddMetaData(object = scrna, metadata = vdj)
  clono.seurat$cdr3 <- as.character(clono.seurat$cdr3)
  clono.seurat$chain <- as.character(clono.seurat$chain)
  clono.seurat$v_gene <- as.character(clono.seurat$v_gene)
  clono.seurat$d_gene <- as.character(clono.seurat$d_gene)
  clono.seurat$j_gene <- as.character(clono.seurat$j_gene)
  
  # Create additional column with only top clones for easier plotting.
  cdr3.count <- as.data.frame(table(as.character(clono.seurat$cdr3)))
  cdr3.count <- cdr3.count[order(-cdr3.count$Freq),]
  rownames(cdr3.count) <- cdr3.count$Var1
  cdr3.count$Var1 <- NULL
  cdr3.count <- cdr3.count[!(rownames(cdr3.count) %in% "None"), , drop = FALSE]
  cdr3.freq <- prop.table(cdr3.count)
  rownames(cdr3.freq) <- rownames(cdr3.count)
  maj <- cdr3.freq[cdr3.freq$Freq > freq.cutoff, , drop = FALSE]
  maj <- rownames(maj)

  clono.seurat$cdr3_group <- "Other"
  clono.seurat$cdr3_group[is.na(clono.seurat$cdr3)] <- "None"
  clono.seurat$cdr3_group[clono.seurat$cdr3 %in% maj] <- 
    clono.seurat$cdr3[clono.seurat$cdr3 %in% maj]
    
  message("Clonotype data added.")
    
  return(clono.seurat)
}
