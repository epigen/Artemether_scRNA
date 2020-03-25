


if(!exists("data.path")) stop("MISSING data.path")
if(!exists("outS")) stop("MISSING outS")
if(!exists("sample.labels")) stop("MISSING sample.labels")
stopifnot(file.exists(data.path))
dir.create(dirout(outS))
(load(data.path))
pbmc.data <- corr.export

if(!exists("seurat.mito.cutoff")) seurat.mito.cutoff <- 0.2
if(!exists("seurat.max.genes")) seurat.max.genes <- Inf
if(!exists("cluster.precision")) cluster.precision <- 0.5



# Apply UMI threshold
cell.umis <- Matrix::colSums(pbmc.data)
if(!is.null(MIN.UMIS)){
  pbmc.data <- pbmc.data[,cell.umis >= MIN.UMIS]
}
cell.genes <- Matrix::colSums(pbmc.data != 0)
if(!is.null(MIN.GENES)){
  pbmc.data <- pbmc.data[,cell.genes >= MIN.GENES]
}

# Set up Seurat Object
pbmc <- CreateSeuratObject(raw.data = pbmc.data, min.cells = 3, min.genes = MIN.GENES, scale.factor=SCALE.FACTOR, project = "X")

# Add labels for aggregated datasets
pbmc@meta.data[["sample"]] <- sample.labels$library_id[as.numeric(gsub("^[A-Z]+\\-(\\d+)$", "\\1", colnames(pbmc@data)))]

# Mitochondrial genes
mito.genes <- grep(pattern = "^MT-", x = rownames(x = pbmc@data), value = TRUE, ignore.case=TRUE)
percent.mito <- Matrix::colSums(pbmc@raw.data[mito.genes, ])/Matrix::colSums(pbmc@raw.data)
pbmc <- AddMetaData(object = pbmc, metadata = percent.mito, col.name = "percent.mito")
ggplot(pbmc@meta.data, aes(y=percent.mito, x=nUMI)) + geom_point(alpha=0.3) + 
  geom_hline(yintercept=seurat.mito.cutoff) + 
  ggtitle(paste(nrow(pbmc@meta.data), "cells"))
ggsave(dirout(outS, "QC_MitoPlot.pdf"))

# Number of genes (many genes probably duplets) ---------------------------
ggplot(pbmc@meta.data, aes(x=nUMI, y=nGene)) + geom_point(alpha=0.3) + 
  geom_hline(yintercept=seurat.max.genes) + ggtitle(paste(nrow(pbmc@meta.data), "cells"))
ggsave(dirout(outS, "QC_GenePlot.pdf"))

# Filter genes
VlnPlot(object = pbmc, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
ggsave(dirout(outS, "QC_Summary.pdf"), height=7, width=7)
pbmc <- FilterCells(object = pbmc, subset.names = c("nGene", "percent.mito"),low.thresholds = c(MIN.GENES, -Inf), high.thresholds = c(seurat.max.genes, seurat.mito.cutoff))
VlnPlot(object = pbmc, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
ggsave(dirout(outS, "QC_Summary2.pdf"), height=7, width=7)

# Normalize and Scale
pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = SCALE.FACTOR)
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

save(pbmc, file=dirout(outS, "SeuratObject",".RData"))


