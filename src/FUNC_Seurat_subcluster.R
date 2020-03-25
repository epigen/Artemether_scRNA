require(methods)
require(Seurat)
require(dplyr)
require(Matrix)
require(pheatmap)
require(data.table)
require(ggplot2)

# REQUIRES
# sample.x - sample name
# data.path - path for input data
# outS - output folder (with dirout)
stopifnot(!is.null(outS))


# OPTIONAL:
#if(!exists("seurat.min.genes")) seurat.min.genes <- 200
#if(!exists("seurat.mito.cutoff")) seurat.mito.cutoff <- 0.15
#if(!exists("seurat.max.genes")) seurat.max.genes <- 5000
# if(!exists("seurat.min.UMIs")) seurat.min.UMIs <- 500
if(!exists("cluster.precision")) cluster.precision <- c(seq(0.5,2.0,0.5))

sessionInfo()

# DEFINE OUTPUT DIRECTORY ------------------------------------------------------------------
dir.create(dirout(outS))
message(" -- saved in -- ", outS)

# ANALYSIS ----------------------------------------------------
# IF RData file exists (preprocessing) just load it
if(!file.exists(dirout(outS, "data.RData"))){
  
  #   stopifnot(!is.null(seurat.mito.cutoff))
  #   stopifnot(!is.null(seurat.max.genes))
  stopifnot(!is.null(cellsToKeep))
  stopifnot(all(cellsToKeep %in% pbmc@cell.names))
  stopifnot(!is.null(pbmc))
  pbmcOrig <- pbmc
  
  
  #   # Apply UMI threshold --> NOT DONE DOESNT WORK, should be  done before in full dataset
  #   cell.umis <- Matrix::colSums(pbmc.data)
  #   if(!is.null(seurat.min.UMIs)){
  #     pbmc.data <- pbmc.data[,cell.umis >= seurat.min.UMIs]
  #   }
  #   cell.genes <- Matrix::colSums(pbmc.data != 0)
  #   if(!is.null(seurat.min.genes)){
  #     pbmc.data <- pbmc.data[,cell.genes >= seurat.min.genes]
  #   }
  
  pbmc <- SubsetData(pbmcOrig, cells.use = cellsToKeep)
  
  # Plot in original tsne ---------------------------------------------------
  tryCatch({
    pDat <- data.table(pbmcOrig@dr$tsne@cell.embeddings,keep.rownames=T)
    pDat$selected <- "no"
    pDat[rn %in% cellsToKeep, selected := "yes"]
    ggplot(pDat, aes(x=tSNE_1, y=tSNE_2, color=selected)) + geom_point(alpha=0.3)
    ggsave(dirout(outS, "SelectedCells.jpg"),height=7, width=7)
  }, error= function(e) message(e))
  
  # Mitochondrial genes
  mito.genes <- grep(pattern = "^MT-", x = rownames(x = pbmc@data), value = TRUE, ignore.case=TRUE)
  percent.mito <- Matrix::colSums(pbmc@raw.data[mito.genes, ])/Matrix::colSums(pbmc@raw.data)
  pbmc <- AddMetaData(object = pbmc, metadata = percent.mito, col.name = "percent.mito")
  ggplot(pbmc@meta.data, aes(y=percent.mito, x=nUMI)) + geom_point(alpha=0.3) + 
    #geom_hline(yintercept=seurat.mito.cutoff) + 
    ggtitle(paste(nrow(pbmc@meta.data), "cells"))
  ggsave(dirout(outS, "QC_MitoPlot.pdf"))
  
  # Number of genes (many genes probably duplets) ---------------------------
  ggplot(pbmc@meta.data, aes(x=nUMI, y=nGene)) + geom_point(alpha=0.3) + 
    #geom_hline(yintercept=seurat.max.genes) + 
    ggtitle(paste(nrow(pbmc@meta.data), "cells"))
  ggsave(dirout(outS, "QC_GenePlot.pdf"))
  
  # Filter genes
  VlnPlot(object = pbmc, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
  ggsave(dirout(outS, "QC_Summary.pdf"), height=7, width=7)
  #pbmc <- FilterCells(object = pbmc, subset.names = c("nGene", "percent.mito"),low.thresholds = c(-Inf, -Inf), high.thresholds = c(Inf, Inf))
  VlnPlot(object = pbmc, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
  ggsave(dirout(outS, "QC_Summary2.pdf"), height=7, width=7)
  
  # Normalize and Scale
  pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
  pbmc <- ScaleData(object = pbmc, vars.to.regress = c("nUMI", "percent.mito"))
  
  # Variable genes
  pbmc <- FindVariableGenes(object = pbmc, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
  #ggsave(dirout(outS, "QC_VariableGenes.pdf"))
  
  # Dimensionality reduction
  pbmc <- RunPCA(pbmc, pc.genes = pbmc@var.genes, do.print = TRUE, pcs.print = 5, genes.print = 5)
  pdf(dirout(outS, "QC_PCs.pdf"), height=20, width=20); PCHeatmap(object = pbmc, pc.use = 1:20, cells.use = 500, do.balanced = TRUE, label.columns = FALSE, use.full = FALSE); dev.off()
  ggplot(data.table(StandDev=pbmc@dr$pca@sdev, PC=1:20), aes(x=PC, y=StandDev)) + geom_point(); ggsave(dirout(outS, "QC_PC_SD.pdf"), height=7, width=7)
  pbmc <- RunTSNE(pbmc, dims.use = 1:10, do.fast = T)
  
  # Clustering
  for(x in cluster.precision){
    pbmc <- FindClusters(pbmc, reduction.type="pca", dims.use = 1:10, resolution = x, print.output = 0, save.SNN = T)
    pbmc <- StashIdent(pbmc, save.name = paste0("ClusterNames_", x))
  }
  
  # WRITE META DATA AND TSNE ------------------------------------------------
  write.table(
    merge(data.table(pbmc@meta.data, keep.rownames=T), data.table(pbmc@dr$tsne@cell.embeddings, keep.rownames=T), by="rn"),
    file=dirout(outS, "MetaData_tSNE.tsv"), sep="\t", quote=F, row.names=F)
  
  save(pbmc, file=dirout(outS,"data.RData"))
  
} else {
  message("Loading processed data")
  
  # LOAD DATA AND UPDATE IF REQUIRED
  load(dirout(outS, "data.RData"))
  update <- FALSE
  
  if(!.hasSlot(pbmc, "version")){
      pbmc <- UpdateSeuratObject(pbmc)
      message("Updating Object Seurat Version")
      update <- TRUE
  }
  
  # Additional clustering if needed
  #   for(x in cluster.precision){
  #     if(is.null(pbmc@meta.data[[paste0("ClusterNames_", x)]])){
  #       update <- TRUE
  #       message("Adding clustering")
  #       pbmc <- FindClusters(pbmc, pc.use = 1:10, resolution = x, print.output = 0, save.SNN = T)
  #       pbmc <- StashIdent(pbmc, save.name = paste0("ClusterNames_", x))
  #     }
  #   }
  
  # SAVE NEW DATASET IF NEEDED
  if(update){
    save(pbmc, file=dirout(outS,"data.RData"))    
  }
}




