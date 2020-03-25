require("project.init")
require(Seurat)
require(dplyr)
require(Matrix)
require(methods)

project.init2("artemether")


R.version

seurat.diff.test <- "negbinom"
seurat.min.genes <- 200

seurat.mito.cutoff <- 0.15
seurat.nGene.cutoff <- 3000

clustering.precision <- seq(0.5, 2.0, 0.5)

out <- "10_Seurat/"
dir.create(dirout(out))

outS <- paste0(out, "DropSeq/")
dir.create(dirout(outS))

cell <- "DropSeq"

if(!file.exists(dirout(outS, cell,".RData"))){
  rawData <- read.csv(dirout(outS, "new.csv"))
  rawData2 <- as.matrix(rawData[,2:ncol(rawData)])
  row.names(rawData2) <- rawData$GENE
  
  pbmc <- CreateSeuratObject(raw.data = rawData2, min.cells = 3, min.genes = seurat.min.genes, 
                             project = "artemether_Dropseq")
  
  mito.genes <- grep(pattern = "^MT-", x = rownames(x = pbmc@data), value = TRUE)
  percent.mito <- Matrix::colSums(pbmc@raw.data[mito.genes, ])/Matrix::colSums(pbmc@raw.data)
  
  pbmc <- AddMetaData(object = pbmc, metadata = percent.mito, col.name = "percent.mito")
  VlnPlot(object = pbmc, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
  ggsave(dirout(outS, "QC_Summary.pdf"), height=7, width=7)
  
  ggplot(pbmc@meta.data, aes(y=percent.mito, x=nUMI)) + geom_point(alpha=0.3) + 
    geom_hline(yintercept=seurat.mito.cutoff) + 
    ggtitle(paste(nrow(pbmc@meta.data), "cells"))
  ggsave(dirout(outS, "QC_MitoPlot.pdf"))
  
  # Number of genes (many genes probably duplets) ---------------------------
  ggplot(pbmc@meta.data, aes(x=nUMI, y=nGene)) + geom_point(alpha=0.3) + 
    geom_hline(yintercept=seurat.nGene.cutoff) + ggtitle(paste(nrow(pbmc@meta.data), "cells"))
  ggsave(dirout(outS, "QC_GenePlot.pdf"))
  
  pbmc <- FilterCells(object = pbmc, subset.names = c("nGene", "percent.mito"),
                      low.thresholds = c(seurat.min.genes, -Inf), high.thresholds = c(seurat.nGene.cutoff, seurat.mito.cutoff))
  
  VlnPlot(object = pbmc, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
  ggsave(dirout(outS, "QC_Summary2.pdf"), height=7, width=7)
  
  pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
  pbmc <- ScaleData(object = pbmc, vars.to.regress = c("nUMI", "percent.mito"))
  
  # PREP DATASET ------------------------------------------------------------
  pbmc <- FindVariableGenes(object = pbmc, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
  ggsave(dirout(outS, "QC_VariableGenes.pdf"))
  pbmc <- RunPCA(pbmc, pc.genes = pbmc@var.genes, do.print = TRUE, pcs.print = 5, genes.print = 5)
  pbmc <- RunTSNE(pbmc, dims.use = 1:10, do.fast = F)
  
  # Clustering
  for(x in clustering.precision){
    pbmc <- FindClusters(pbmc, reduction.type="pca", dims.use = 1:10, resolution = x, print.output = 0, save.SNN = T)
    pbmc <- StashIdent(pbmc, save.name = paste0("ClusterNames_", x))
  }
  
  pbmc@meta.data$sample <- gsub("(.+?)\\.(N|A|C|T|G)+", "\\1", colnames(pbmc@data))
  
  save(pbmc, file=dirout(outS, cell,".RData"))
} else {
  load(file=dirout(outS, cell,".RData"))
  
  update <- FALSE
  for(x in clustering.precision){
    if(is.null(pbmc@meta.data[[paste0("ClusterNames_", x)]])){
      update <- TRUE
      pbmc <- FindClusters(pbmc, reduction.type="pca", dims.use = 1:10, resolution = x, print.output = 0, save.SNN = T)
      pbmc <- StashIdent(pbmc, save.name = paste0("ClusterNames_", x))
    }
  }
  
  if(update){
    save(pbmc, file=dirout(outS, cell,".RData"))
  }
}



max.nr.of.cores <- 4
extra.genes.to.plot <- unique(fread("metadata//PancreasMarkers.tsv")$Human_GeneName)
source(paste0(Sys.getenv("CODEBASE"), "10x_datasets/src/FUNC_Enrichr.R"))
source(paste0(Sys.getenv("CODEBASE"), "10x_datasets/src/FUNC_Seurat2.R"), echo=TRUE)
