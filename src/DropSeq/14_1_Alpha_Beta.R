require("project.init")
require(Seurat)
require(Matrix)
require(methods)
project.init2("artemether")

inDir <- 
outS <- "14_1_Alpha_Beta/"
dir.create(dirout(outS))

cluster.precision <- c(0.5, 1.0, 1.5, 2.0)

cell <- "Alpha_Beta"

if(!file.exists(dirout(outS, cell,".RData"))){
  (load(file=dirout("10_Seurat/", "DropSeq/","DropSeq.RData")))
  pDat <- fread(dirout("13_2_CellTypes_noGABA/", "MetaData.tsv"))
  
  stopifnot(all(row.names(subset(pbmc@meta.data, sample %in% c("Artemether", "DMSO"))) == pDat$cell))
  
  pbmc@meta.data$Celltype <- pDat[match(pbmc@cell.names, pDat$cell)]$Celltype
  
  pbmcOrig <- pbmc
  
  print(unique(pbmcOrig@meta.data$Celltype))
  
  # SUBSET DATA -------------------------------------------------------------
  str(cellToKeep <- row.names(subset(pbmc@meta.data, Celltype %in% c("INS+", "GCG+", "INS+ GCG+"))))
  
  dir.create(dirout(outS))
  
  pbmc <- SubsetData(pbmcOrig, cells.use = cellToKeep)
  
  pDat <- data.table(pbmcOrig@dr$tsne@cell.embeddings,keep.rownames=T)
  pDat$selected <- "no"
  pDat[rn %in% cellToKeep, selected := "yes"]
  ggplot(pDat, aes(x=tSNE_1, y=tSNE_2, color=selected)) + geom_point(alpha=0.3)
  ggsave(dirout(outS, "SelectedCells.pdf"))
  
  # PREP DATASET ------------------------------------------------------------
  pbmc <- FindVariableGenes(pbmc, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
  pbmc <- RunPCA(pbmc, pc.genes = pbmc@var.genes, do.print = TRUE, pcs.print = 5, genes.print = 5)
  pbmc <- RunTSNE(pbmc, dims.use = 1:10, do.fast = F)
  
  # Clustering
  for(x in cluster.precision){
    pbmc <- FindClusters(pbmc, reduction.type="pca", dims.use = 1:10, resolution = x, print.output = 0, save.SNN = T)
    pbmc <- StashIdent(pbmc, save.name = paste0("ClusterNames_", x))
  }
  pbmc@meta.data[["Alpha_Beta"]] <- NULL
  save(pbmc, file=dirout(outS, cell,".RData"))
} else {
  load(file=dirout(outS, cell,".RData"))
  
  update <- FALSE
  for(x in cluster.precision){
    if(is.null(pbmc@meta.data[[paste0("ClusterNames_", x)]])){
      update <- TRUE
      pbmc <- FindClusters(pbmc, reduction.type="pca", dims.use = 1:10, resolution = x, print.output = 0, save.SNN = T)
      pbmc <- StashIdent(pbmc, save.name = paste0("ClusterNames_", x))
    }
  }
  
  resNames <- colnames(pbmc@meta.data)[grepl("res\\.", colnames(pbmc@meta.data))]
  for(rr in resNames[! resNames %in% paste0("res.", cluster.precision)]){
    pbmc@meta.data[[rr]] <- NULL
    update <- TRUE
  }
  if(update){
    save(pbmc, file=dirout(outS, cell,".RData"))
  }
}

pDat <- data.table(pbmc@meta.data, pbmc@dr$pca@cell.embeddings, SST = pbmc@data["SST",])
ggplot(pDat, aes(x=PC1, y=PC2, color=Celltype)) + geom_point()
ggplot(pDat, aes(x=PC1, y=PC2, color=SST)) + geom_point(alpha=0.5) +scale_color_gradient(low="lightgrey", high="blue")


extra.genes.to.plot <- unique(fread("metadata//PancreasMarkers.tsv")$Human_GeneName)
source(paste0(Sys.getenv("CODEBASE"), "10x_datasets/src/FUNC_Seurat2.R"), echo=TRUE)

(load(paste(Sys.getenv("CODEBASE"), "slice/data/hs_km.Rda", sep="")))
SLICE.km <- km
SLICE.cellIdentity <- factor(pbmc@meta.data$sample)
source(paste0(Sys.getenv("CODEBASE"), "slice/slice.R"))
source(paste0(Sys.getenv("CODEBASE"), "10x_datasets/src/FUNC_SLICE2.R"), echo=TRUE)

Monocle.metadata.fields <- c("Slice2")
source(paste0(Sys.getenv("CODEBASE"), "10x_datasets/src/FUNC_Monocle.R"), echo=TRUE)

extra.genes.to.plot <- unique(fread("metadata//PancreasMarkers.tsv")$Human_GeneName)
source(paste0(Sys.getenv("CODEBASE"), "10x_datasets/src/FUNC_Seurat2.R"), echo=TRUE)
