require("project.init")
require(Seurat)
require(Matrix)
require(methods)
project.init2("artemether")

source("src/FUNC_Enrichr.R"))

inDir <- paste0("10_H_01_Seurat/", "Subcluster/")

load(dirout(inDir, "data",".RData"))
pbmcOLD <- pbmc

out <- "30_01_ClusteringExperiments/"
dir.create(dirout(out))

# SUBCLUSTER --------------------------------------------------------------
meta <- loadHMeta2()
meta <- meta[rn %in% pbmc@cell.names]
cellsToKeep.list <- with(meta, split(rn, factor(paste(treatment, replicate, sep="_"))))
str(cellsToKeep.list)

xnam <- names(cellsToKeep.list)[1]
for(xnam in names(cellsToKeep.list)){
  # do this for each
  outS <- paste0(out, xnam, "/")
  cellsToKeep <- cellsToKeep.list[[xnam]]
  cluster.precision <- 0.5
  sample.x <- xnam
  pbmc <- pbmcOLD
  source("src/FUNC_Seurat_subcluster.R",echo=TRUE)
  ggplot(data.table(pbmc@dr$tsne@cell.embeddings), aes(x=tSNE_1, y=tSNE_2)) + geom_point(alpha=0.1)
  ggsave(dirout(outS, "TSNE.jpg"), w=4,h=4)
}

